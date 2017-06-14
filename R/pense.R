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
#'      See \code{\link{en_options}} for details.
#' @param initest_options additional options for computing the initial
#'      estimator. Ignored if \code{initial = "warm0"}.
#'      See \code{\link{initest_options}} for details.
#'
#' @return An object of class \code{"pense"} with elements
#' \item{call}{the call that produced this object.}
#' \item{lambda_opt}{the optimal value of the regularization parameter
#'      according to CV.}
#' \item{coefficients}{the vector of coefficients for the optimal lambda
#'      \code{lambda_opt}.}
#' \item{residuals}{the residuals, that is response minus fitted values for
#'      \code{lambda_opt}.}
#' \item{scale}{the estimated scale for \code{lambda_opt}.}
#' \item{lambda_grid}{a matrix with values of lambda in the first
#'      column, the \code{cv_objective} as second column and serveral
#'      statistics of the solution in the following columns.}
#' \item{alpha,standardize,pense_options,en_options,initest_options}{
#'      the given arguments.}
#' @export
#'
#' @importFrom stats mad median
#' @importFrom robustbase scaleTau2
pense <- function(X, y,
                  alpha = 0.5,
                  nlambda = 100, lambda = NULL, lambda_min_ratio = NULL,
                  standardize = TRUE,
                  initial = c("warm", "warm0", "cold"),
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

    scale_x <- 1
    mux <- 0
    muy <- 0

    ## Standardize data -- only needed to compute the grid of lambda values
    if (isTRUE(standardize)) {
        scale_x <- apply(X, 2, mad)
        mux <- apply(X, 2, median)
        muy <- median(y)

        Xs <- scale(X, center = mux, scale = scale_x)
        yc <- y - muy

        if (!is.null(lambda)) {
            lambda <- lambda / max(scale_x)
        }
    } else {
        Xs <- X
        yc <- y
    }

    ## Generate grid of lambda-values
    if (is.null(lambda)) {
        lambda <- build_lambda_grid(
            Xs,
            yc,
            alpha,
            nlambda,
            standardize = FALSE,
            lambda_min_ratio = lambda_min_ratio
        )
    }

    ## Reverse the order if we use a warm-0 start
    lambda <- sort(lambda, decreasing = isTRUE(initial == "warm0"))

    ## CV function (needs X, y, and scale_x available in the environment)
    pense_job_cv <- function(job, ...) {
        if (length(job$segment) == 0L) {
            X_train <- X
            y_train <- y
        } else {
            X_train <- X[-job$segment, , drop = FALSE]
            y_train <- y[-job$segment]

        }
        est.all <- pense_coldwarm(
            X = X_train,
            y = y_train,
            lambda_grid = job$lambda,
            ...
        )

        residuals <- NULL

        if (length(job$segment) == 0L) {
            residuals <- vapply(est.all, function(est) {
                est$resid
            }, FUN.VALUE = y, USE.NAMES = FALSE)
        } else {
            X_test <- X[job$segment, , drop = FALSE]
            y_test <- y[job$segment]
            residuals <- vapply(est.all, function(est) {
                y_test - est$intercept - X_test %*% est$beta
            }, FUN.VALUE = numeric(length(job$segment)), USE.NAMES = FALSE)
        }

        sol_stats <- vapply(est.all, function (est) {
            c(
                objF = est$objF,
                scale = est$scale,
                beta_L1 = sum(abs(est$beta / scale_x)),
                beta_L2 = sqrt(sum((est$beta / scale_x)^2))
            )
        }, FUN.VALUE = numeric(4L), USE.NAMES = TRUE)

        return(list(
            residuals = residuals,
            sol_stats = sol_stats
        ))
    }

    ## Perform CV (if we have more than a single lambda value)
    if(nlambda > 1L) {
        cv_objective_fun <- if (missing(cv_objective)) {
            scaleTau2
        } else {
            match.fun(cv_objective)
        }

        # Create CV segments
        cv_segments <- if (cv_k > 1L) {
            split(seq_len(dX[1L]), sample(rep_len(seq_len(cv_k), dX[1L])))
        } else {
            list(integer(0L))
        }

        # Define CV jobs (i.e., all combinations of CV splits and warm resets)
        subgrid_lengths <- nlambda %/% warm_reset
        overlength <- nlambda %% warm_reset
        subgrid_lengths <- c(rep.int(subgrid_lengths, warm_reset - overlength),
                             rep.int(subgrid_lengths + 1L, overlength))

        lambda_subgrids <- split(lambda,
                                 rep.int(seq_len(warm_reset), subgrid_lengths))

        jobgrid <- lapply(cv_segments, function(cvs) {
            lapply(lambda_subgrids, function(lvs) {
                list(
                    segment = cvs,
                    lambda = lvs
                )
            })
        })

        jobgrid <- unlist(jobgrid, recursive = FALSE, use.names = FALSE)

        # Setup cluster
        cluster <- setupCluster(
            ncores,
            cl,
            export = c("X", "y", "scale_x"),
            eval = {
                library(pense)
            }
        )

        # Run all jobs (combination of all CV segments and all lambda-grids)
        tryCatch({
            cv_results <- cluster$lapply(
                jobgrid,
                pense_job_cv,
                alpha = alpha,
                standardize = standardize,
                start_0 = isTRUE(initial == "warm0"),
                pense_options = options,
                initest_options = init_options,
                en_options = en_options
            )

            cv_results <- split(cv_results, rep.int(seq_len(warm_reset), cv_k))
        },
        finally = {
            cluster$stopCluster()
        })

        # Collect all prediction errors for each lambda sub-grid and determine the optimal lambda
        cv_performance <- unlist(lapply(cv_results, function(cv_res) {
            pred_resids <- do.call(rbind, lapply(cv_res, "[[", "residuals"))
            apply(pred_resids, 2, cv_objective_fun)
        }), use.names = FALSE)

        cv_stats <- do.call(rbind, lapply(cv_results, function (cv_res) {
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
        }))

        lambda_grid <- data.frame(
            lambda = lambda * max(scale_x),
            cv_performance = cv_performance,
            cv_stats
        )
        lambda_opt <- lambda[which.min(cv_performance)]
    } else {
        lambda_grid <- NULL
        lambda_opt <- lambda[1L]
    }

    ## Fit PENSE with optimal lambda using a cold start
    opt_est <- pense_coldwarm(
        X = Xs,
        y = yc,
        alpha = alpha,
        lambda_grid = lambda_opt,
        standardize = FALSE,
        start_0 = isTRUE(initial == "warm0"),
        pense_options = options,
        initest_options = init_options,
        en_options = en_options
    )[[1L]]

    ## Un-standardize the coefficients
    if (isTRUE(standardize)) {
        opt_est$beta <- opt_est$beta / scale_x
        opt_est$intercept <- opt_est$intercept + muy - drop(opt_est$beta %*% mux)
    }

    return(structure(list(
        residuals = opt_est$resid,
        coefficients = nameCoefVec(c(opt_est$intercept, opt_est$beta), X),
        objective = opt_est$objF,
        lambda_opt = lambda_opt * max(scale_x),
        lambda_grid = lambda_grid,
        scale = opt_est$scale,
        alpha = alpha,
        standardize = standardize,
        pense_options = options,
        initest_options = init_options,
        en_options = en_options,
        call = call
    ), class = "pense"))
}

