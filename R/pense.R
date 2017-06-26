#' Penalized Elasitc Net S-estimators for Regression
#'
#' @section Parallelization:
#' With the parameter \code{ncores}, the number of available processor cores
#' can be set. The grid of lambda values is split into \code{ncores} (almost)
#' equally sized chunks. On each core, the cold initial estimate is computed
#' for the lower endpoint of the lambda sequence. The fully iterated PENSE
#' estimate is then used as warm-start for the next lambda value and so on.
#' Therefore, if \code{ncores = 1}, only one cold initial estimate is computed
#' and all subsequent values of lambda take the previous estimate as
#' warm-start. Similarly, if \code{ncores = nlambda}, a cold initial estimate
#' is computed for every lambda on the grid and no warm-starts are necessary.
#'
#' @param X design matrix with predictors.
#' @param y response vector.
#' @param alpha elastic net mixing parameter with \eqn{0 \le \alpha \le 1}.
#'      \code{alpha = 1} is the LASSO penatly, and \code{alpha = 0} the Ridge
#'      penalty.
#' @param nlambda if \code{lambda} is not given or \code{NULL} (default),
#'      a grid of \code{nlambda} lambda values is generated based on the data.
#' @param lambda a single value or a grid of values for the regularization
#'      parameter lambda.
#'      Assumed to be on the same scale as the data and adjusted for
#'      S-estimation. Defaults to \code{NULL}, which means a grid of lambda
#'      values is automatically generated (see parameter \code{nlambda}).
#'      If given and \code{standardize = TRUE}, the lambda values will be
#'      adjusted accordingly.
#' @param lambda_min_ratio If the grid should be chosen automatically, the
#'      ratio of the smallest lambda to the (computed) largest lambda.
#' @param standardize should the data be standardized robustly? Estimates
#'      are returned on the original scale. Defaults to \code{TRUE}.
#' @param cv_k number of cross-validation segements to use to choose the optimal
#'      lambda from the grid. If only a single value of lambda is given,
#'      cross-validation can still done to estimate the prediction performance
#'      at this particular lambda.
#' @param cv_objective a function (name) to compute the CV performance.
#'      By default, the robust tau-scale is used.
#' @param initial how to initialize the estimator at a new lambda in the grid.
#'      The default, \code{"warm"}, computes a cold initial estimator at several
#'      lambda values and uses the PENSE coefficient to warm-start the
#'      estimator at the next larger lambda value. A variant, \code{"warm0"},
#'      will initialze PENSE at the largest lambda value with an all-0
#'      coefficient vector and use the PENSE result for the next smaller
#'      lambda value. \code{"cold"} computes the full initial estimator at
#'      every lambda value.
#' @param warm_reset if \code{initial = "warm"}, how often should the warm-start
#'      be reset to a cold initial estimate?
#' @param ncores,cl the number of processor cores or an actual parallel cluster
#'      to use to estimate the optimal value of lambda. See details for more
#'      information.
#' @param options additional options for the PENSE algorithm.
#'      See \code{\link{pense_options}} for details.
#' @param en_options additional options for the EN algorithm.
#'      See \code{\link{elnet}} and \code{\link{en_options}} for details.
#' @param initest_options additional options for computing the initial
#'      estimator. Ignored if \code{initial = "warm0"}.
#'      See \code{\link{initest_options}} for details.
#'
#' @return An object of class \code{"pense"} with elements
#' \item{call}{the call that produced this object.}
#' \item{lambda_opt}{the optimal value of the regularization parameter
#'      according to CV.}
#' \item{coefficients}{a matrix of coefficients for each lambda in the grid.}
#' \item{residuals}{a matrix of residuals for each lambda in the grid.}
#' \item{scale}{the estimated scales each lambda}
#' \item{cv_lambda_grid}{a data frame with values of lambda in the first
#'      column, and the CV error as well as serveral
#'      statistics of the solution in the following columns.}
#' \item{alpha,standardize,pense_options,en_options,initest_options}{
#'      the given arguments.}
#' @export
#'
#' @importFrom stats mad median
#' @importFrom robustbase scaleTau2
#' @importFrom Matrix norm drop
pense <- function(X, y,
                  alpha = 0.5,
                  nlambda = 50, lambda = NULL, lambda_min_ratio = NULL,
                  standardize = TRUE,
                  initial = c("warm0", "warm", "cold"),
                  warm_reset = 10,
                  cv_k = 5, cv_objective,
                  ncores = getOption("mc.cores", 1L), cl = NULL,
                  options = pense_options(),
                  init_options = initest_options(),
                  en_options = en_options_aug_lars()) {
    ##
    ## First all arguments are being sanity-checked
    ##
    y <- drop(y)

    dY <- dim(y)
    yl <- length(y)
    dX <- dim(X)

    if (is.data.frame(X)) {
        X <- data.matrix(X)
    }

    if (!is.matrix(X) || !is.numeric(X)) {
        stop("`X` must be a numeric matrix")
    }

    if (is.null(yl) || (!is.null(dY) && length(dY) != 1L) || !is.numeric(y)) {
        stop("`yl` must be a numeric vector")
    }

    if (dX[1L] != yl) {
        stop("The number of observations in `X` and `y` does not match")
    }

    if (any(!is.finite(X))) {
        stop("`X` must not contain infinite, NA, or NaN values")
    }

    if (any(!is.finite(y))) {
        stop("`y` must not contain infinite, NA or NaN values")
    }

    alpha <- .check_arg(alpha, "numeric", range = c(0, 1),
                        range_test_lower = ">=", range_test_upper = "<=")

    standardize <- .check_arg(standardize, "logical")
    cv_k <- .check_arg(cv_k, "integer", range = 0)

    initial <- match.arg(initial)

    warm_reset <- if (isTRUE(initial == "warm")) {
        if (is.null(warm_reset)) {
            warm_reset <- nlambda
        }

        .check_arg(warm_reset, "integer", range = 0)
    } else if (isTRUE(initial == "cold")) {
        Inf
    } else {
        1L
    }

    if (!is.null(lambda)) {
        lambda <- .check_arg(lambda, "numeric", range = 0, length = NULL)
        nlambda <- length(lambda)
    } else {
        nlambda <- .check_arg(nlambda, "integer", range = 0)
    }
    warm_reset <- min(nlambda, warm_reset)

    if (is.null(lambda) && !is.null(lambda_min_ratio)) {
        lambda_min_ratio <- .check_arg(lambda_min_ratio, "numeric",
                                       range = c(0, 1))
    }

    ##
    ## Sanity checks for some arguments. Some combinations are known to
    ## be performing poorly.
    ##
    if (initial != "warm0" && alpha < 1 && init_options$pscMethod == "rr" &&
        en_options$algorithm == 1L) {
        message("Approximation of PSCs does not work well with the DAL ",
                "algorithm for EN due to the large number of observations ",
                "in the augmented data matrix. ",
                'Consider using `psc_method = "exact"` or use the augmented ',
                "LARS algorithm for EN.")
    }

    if (en_options$algorithm == 1L && dX[1L] > 500L && dX[2L] < 1000L) {
        message("The DAL algorithm for elastic net scales with the number of ",
                "observations. Given the number of variables, augmented LARS ",
                "might be the faster option.")
    }

    ## store the call
    call <- match.call()
    call[[1L]] <- as.name("pense")

    std_data <- standardize_data(X, y, standardize)
    Xs <- std_data$xs
    yc <- std_data$yc
    scale_x <- std_data$scale_x

    ## Generate grid of lambda-values
    if (is.null(lambda)) {
        lambda <- build_lambda_grid(
            Xs,
            yc,
            alpha,
            nlambda,
            lambda_min_ratio = lambda_min_ratio
        )
    }

    lambda <- lambda / max(scale_x)

    ## Reverse the order if we use a warm-0 start
    lambda <- sort(lambda, decreasing = isTRUE(initial == "warm0"))

    ## Setup cluster
    cluster <- setupCluster(
        ncores,
        cl,
        export = c("X", "y", "scale_x"),
        eval = {
            library(pense)
        }
    )

    ## Function returning the estimates and residuals at every lambda value
    ## (needs X, y, and scale_x available in the environment)
    pense_est_job <- function(job, ...) {
        if (length(job$segment) == 0L) {
            X_train <- X
            y_train <- y
            ret_coefs <- TRUE
        } else {
            X_train <- X[-job$segment, , drop = FALSE]
            y_train <- y[-job$segment]
            ret_coefs <- FALSE
        }

        est_all <- pense_coldwarm(
            X = X_train,
            y = y_train,
            lambda_grid = job$lambda,
            ...
        )

        residuals <- NULL
        intercept <- NULL
        beta <- NULL
        if (isTRUE(ret_coefs)) {
            beta <- do.call(cbind, lapply(est_all, "[[", "beta"))
            intercept <- unlist(lapply(est_all, "[[", "intercept"))
        }

        if (length(job$segment) == 0L) {
            residuals <- vapply(
                est_all,
                function(est) {
                    est$resid
                },
                FUN.VALUE = y, USE.NAMES = FALSE
            )
        } else {
            X_test <- X[job$segment, , drop = FALSE]
            y_test <- y[job$segment]
            residuals <- vapply(
                est_all,
                function(est) {
                    drop(y_test - est$intercept - X_test %*% est$beta)
                },
                FUN.VALUE = numeric(length(job$segment)), USE.NAMES = FALSE)
        }

        sol_stats <- vapply(
            est_all,
            function (est) {
                c(
                    objF = est$objF,
                    scale = est$scale,
                    beta_L1 = norm(est$beta, "1"),
                    beta_L2 = norm(est$beta, "F")
                )
            },
            FUN.VALUE = numeric(4L), USE.NAMES = TRUE)

        return(list(
            residuals = residuals,
            intercept = intercept,
            beta = beta,
            sol_stats = sol_stats
        ))
    }

    ##
    ## Define the lambda sub-grids
    ##
    lambda_subgrid_lengths <- nlambda %/% warm_reset
    overlength <- nlambda %% warm_reset
    lambda_subgrid_lengths <- c(
        rep.int(lambda_subgrid_lengths, warm_reset - overlength),
        rep.int(lambda_subgrid_lengths + 1L, overlength)
    )

    lambda_subgrids <- split(
        lambda,
        rep.int(seq_len(warm_reset), lambda_subgrid_lengths)
    )

    ## The jobs for the "full" (original) data set
    jobgrid_full <- lapply(lambda_subgrids, function(lvs) {
        list(
            segment = integer(0L),
            lambda = lvs
        )
    })

    ## Perform CV (if we have more than a single lambda value)
    if((nlambda > 1L) && (cv_k > 1L))  {
        # Create CV segments
        cv_segments <- split(
            seq_len(dX[1L]),
            sample(rep_len(seq_len(cv_k), dX[1L]))
        )

        # Define CV jobs (i.e., all combinations of CV splits and
        # lambda subgrids)
        jobgrid_cv <- lapply(cv_segments, function(cvs) {
            lapply(lambda_subgrids, function(lvs) {
                list(
                    segment = cvs,
                    lambda = lvs
                )
            })
        })

        jobgrid_cv <- unlist(jobgrid_cv, recursive = FALSE, use.names = FALSE)
    } else {
        jobgrid_cv <- list()
    }

    # Run all jobs (original & CV)
    all_jobs <- c(jobgrid_full, jobgrid_cv)
    tryCatch({
        all_results <- cluster$lapply(
            all_jobs,
            pense_est_job,
            alpha = alpha,
            standardize = standardize,
            start_0 = isTRUE(initial == "warm0"),
            pense_options = options,
            initest_options = init_options,
            en_options = en_options
        )

        full_results <- all_results[seq_along(jobgrid_full)]
        cv_results <- all_results[-seq_along(jobgrid_full)]
    },
    finally = {
        cluster$stopCluster()
    })

    ## Gather all coefficients and residuals from the full results
    int_ests <- unlist(lapply(full_results, "[[", "intercept"))
    beta_ests <- do.call(cbind, lapply(full_results, "[[", "beta"))
    residuals <- unlist(
        lapply(full_results, "[[", "residuals"),
        use.names = FALSE
    )
    residuals <- matrix(residuals, nrow = dX[1L], ncol = nlambda)
    scale_ests <- unlist(
        lapply(full_results, function (lambda_res) {
            lambda_res$sol_stats["scale", ]
        }),
        use.names = FALSE
    )
    obj_fun_vals <- unlist(
        lapply(full_results, function (lambda_res) {
            lambda_res$sol_stats["objF", ]
        }),
        use.names = FALSE
    )

    ## Gather all stats and results from the CV jobs
    if (length(cv_results) > 0L) {
        cv_results <- split(cv_results, rep.int(seq_len(warm_reset), cv_k))

        cv_objective_fun <- if (missing(cv_objective)) {
            scaleTau2
        } else {
            match.fun(cv_objective)
        }

        # Collect all prediction errors for each lambda sub-grid and determine the optimal lambda
        cv_obj <- t(do.call(cbind, lapply(
            cv_results,
            function(cv_res) {
                cv_obj_fvals <- unlist(lapply(cv_res, function (lambda_res) {
                    apply(lambda_res$residuals, 2, cv_objective_fun)
                }))
                cv_obj_fvals <- matrix(cv_obj_fvals, ncol = length(cv_res))
                apply(cv_obj_fvals, 1, function (x) {
                    c(cvavg = mean(x), cvsd = sd(x))
                })

                # pred_resids <- do.call(rbind, lapply(cv_res, "[[", "residuals"))
                # apply(pred_resids, 2, cv_objective_fun)
            }))
        )

        cv_stats <- do.call(rbind, lapply(
            cv_results,
            function (cv_res) {
                sol_stats <- unlist(lapply(cv_res, "[[", "sol_stats"))
                sol_stats <- array(
                    sol_stats,
                    dim = c(4L, length(sol_stats) %/% (cv_k * 4L), cv_k),
                    dimnames = list(
                        c("obj_fun", "s_scale", "beta_L1", "beta_L2"),
                        NULL,
                        NULL
                    )
                )
                apply(sol_stats, 1L, rowMeans, na.rm = TRUE)
            })
        )

        cv_lambda_grid <- data.frame(
            lambda = lambda * max(scale_x),
            cv_obj,
            cv_stats
        )
        lambda_opt <- lambda[which.min(cv_obj[, "cvavg"])]
    } else {
        cv_lambda_grid <- NULL
        lambda_opt <- lambda[1L]
    }

    ## Order results on lambda grid by increasing lambda
    lambda_order <- sort.list(lambda, method = "quick", na.last = NA)
    lambda <- lambda[lambda_order]
    coef_ests <- rbind(
        int_ests[lambda_order],
        beta_ests[ , lambda_order, drop = FALSE]
    )
    residuals <- residuals[, lambda_order, drop = FALSE]
    scale_ests <- scale_ests[lambda_order]
    obj_fun_vals <- obj_fun_vals[lambda_order]
    cv_lambda_grid <- cv_lambda_grid[lambda_order, , drop = FALSE]

    return(structure(list(
        residuals = residuals,
        coefficients = nameCoefVec(coef_ests, X),
        lambda = lambda * max(scale_x),
        scale = scale_ests,
        objective = obj_fun_vals,
        cv_lambda_grid = cv_lambda_grid,
        lambda_opt = lambda_opt * max(scale_x),
        alpha = alpha,
        standardize = standardize,
        pense_options = options,
        initest_options = init_options,
        en_options = en_options,
        call = call
    ), class = "pense"))
}
