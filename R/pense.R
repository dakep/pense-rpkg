#' Penalized Elastic Net S-estimators for Regression
#'
#' Computes the highly robust Penalized Elastic Net S-estimators (PENSE) for
#' linear regression models.
#'
#' The PENSE estimate minimizes the robust M-scale of the residuals penalized
#' by the L1 and L2 norm of the regression coefficients (elastic net penalty).
#' The level of penalization is chosen to minimize the \code{cv_k}-fold
#' cross-validated prediction error (using a robust measure).
#'
#' @section Initial Estimate:
#' By default (\code{initial == "warm"}), the method does not compute a
#' full initial estimate at each
#' lambda value in the grid, but only at \code{warm_reset} of the lambda
#' values. At the remaining lambda values, the estimate at the previous
#' lambda value is used to initialize the estimator (the lambda grid is
#' first traversed in descending and then in ascending direction). If
#' \code{warm_reset} is 1, only the 0-vector is used to initialize PENSE at the
#' largest penalty value. No further initial estimates are computed.
#'
#' If \code{initial == "cold"}, a full initial estimate is computed at each
#' lambda value. This is equal to setting \code{warm_reset} to
#' \code{length(lambda)}.
#'
#' @param x design matrix with predictors.
#' @param y response vector.
#' @param alpha elastic net mixing parameter with \eqn{0 \le \alpha \le 1}.
#'      \code{alpha = 1} is the LASSO penalty, and \code{alpha = 0} the Ridge
#'      penalty.
#' @param nlambda if \code{lambda} is not given or \code{NULL} (default),
#'      a grid of \code{nlambda} lambda values is generated based on the data.
#' @param lambda a single value or a grid of values for the regularization
#'      parameter lambda.
#'      Assumed to be on the same scale as the data and adjusted for
#'      S-estimation. If missing a grid of lambda
#'      values is automatically generated (see parameter \code{nlambda}).
#'      If given and \code{standardize = TRUE}, the lambda values will be
#'      adjusted accordingly.
#' @param lambda_min_ratio If the grid should be chosen automatically, the
#'      ratio of the smallest lambda to the (computed) largest lambda.
#' @param standardize should the data be standardized robustly? Estimates
#'      are returned on the original scale. Defaults to \code{TRUE}.
#' @param cv_k number of cross-validation segments to use to choose the optimal
#'      lambda from the grid. If only a single value of lambda is given,
#'      cross-validation can still done to estimate the prediction performance
#'      at this particular lambda.
#' @param cv_objective a function (name) to compute the CV performance.
#'      By default, the robust tau-scale is used.
#' @param initial how to initialize the estimator at a new lambda in the grid.
#'      The default, \code{"warm"}, computes a cold initial estimator at several
#'      lambda values and uses the PENSE coefficient to warm-start the
#'      estimator at the next larger lambda value. At the largest value in
#'      the lambda grid, PENSE will be initialized with the 0-vector.
#'      \code{"cold"} computes the full initial estimator at
#'      every lambda value.
#' @param warm_reset if \code{initial = "warm"} (default), how many cold initial
#'      estimates be computed?
#' @param ncores,cl the number of processor cores or an actual parallel cluster
#'      to use to estimate the optimal value of lambda. See
#'      \code{\link{makeCluster}} on how to create a cluster manually.
#' @param options additional options for the PENSE algorithm.
#'      See \code{\link{pense_options}} for details.
#' @param en_options additional options for the EN algorithm.
#'      See \code{\link{elnet}} and \code{\link{en_options}} for details.
#' @param init_options additional options for computing the cold initial
#'      estimates.
#'      Ignored if \code{initial = "warm"} and \code{warm_reset = 0}.
#'      See \code{\link{initest_options}} for details.
#'
#' @return An object of class \code{"pense"} with elements
#' \item{lambda}{grid of regularization parameter values for which an estimate
#'      is available.}
#' \item{lambda_opt}{the optimal value of the regularization parameter
#'      according to CV.}
#' \item{coefficients}{a sparse matrix of coefficients for each lambda in the
#'      grid.}
#' \item{residuals}{a matrix of residuals for each lambda in the grid.}
#' \item{cv_lambda_grid}{a data frame with CV prediction errors and several
#'      statistics of the solutions.}
#' \item{scale}{the estimated scales each lambda in the grid.}
#' \item{objective}{value of the objective function at each lambda in the grid.}
#' \item{adjusted}{necessary information to compute the corrected EN estimates.}
#' \item{call}{the call that produced this object.}
#' \item{...}{values of the given arguments.}
#'
#' @seealso To improve the S-estimate with an M-step, see \code{\link{pensem}}.
#'
#' @example examples/pense-1.R
#'
#' @export
#'
#' @importFrom stats mad median weighted.mean
#' @importFrom robustbase scaleTau2
#' @importFrom Matrix norm drop Diagonal colSums
#' @useDynLib pense, .registration = TRUE
pense <- function(x, y,
                  alpha = 0.5,
                  nlambda = 50, lambda, lambda_min_ratio,
                  standardize = TRUE,
                  initial = c("warm", "cold"),
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
    dx <- dim(x)

    if (is.data.frame(x)) {
        x <- data.matrix(x)
    }

    if (!is.matrix(x) || !is.numeric(x)) {
        stop("`x` must be a numeric matrix")
    }

    if (is.null(yl) || (!is.null(dY) && length(dY) != 1L) || !is.numeric(y)) {
        stop("`yl` must be a numeric vector")
    }

    if (dx[1L] != yl) {
        stop("The number of observations in `x` and `y` does not match")
    }

    if (any(!is.finite(x))) {
        stop("`x` must not contain infinite, NA, or NaN values")
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
        .check_arg(warm_reset, "integer", range = 0)
    } else if (isTRUE(initial == "cold")) {
        Inf
    } else {
        1L
    }

    if (!missing(lambda) && !is.null(lambda)) {
        lambda <- .check_arg(lambda, "numeric", range = 0, length = NULL)
        nlambda <- length(lambda)
    } else {
        nlambda <- .check_arg(nlambda, "integer", range = 0)
        lambda <- NULL
    }
    warm_reset <- min(nlambda, warm_reset)

    if (missing(lambda_min_ratio) || is.null(lambda_min_ratio)) {
        lambda_min_ratio <- .default_lambda_min_ratio(x)
    }

    if (is.null(lambda) && !is.null(lambda_min_ratio)) {
        lambda_min_ratio <- .check_arg(
            lambda_min_ratio,
            "numeric",
            range = c(0, 1)
        )
    }

    ##
    ## Sanity checks for some arguments. Some combinations are known to
    ## be performing poorly.
    ##
    if (!(initial == "warm" && warm_reset == 1L) &&
        alpha > 0 && alpha < 1 &&
        init_options$pscMethod == "rr" &&
        en_options$algorithm == 1L) {
        message("Approximation of PSCs does not work well with the DAL ",
                "algorithm for EN due to the large number of observations ",
                "in the augmented data matrix. ",
                'Consider using `initest_options(psc_method = "exact")` or ',
                "use the augmented LARS algorithm for EN.")
    }

    if (en_options$algorithm == 1L && dx[1L] > 500L && dx[2L] < 1000L) {
        message("The DAL algorithm for elastic net scales with the number of ",
                "observations. Given the number of variables, augmented LARS ",
                "might be the faster option.")
    }

    ## store the call
    call <- match.call()
    call[[1L]] <- as.name("pense")

    ## Generate grid of lambda-values
    if (is.null(lambda)) {
        lambda <- build_lambda_grid(
            x,
            y,
            alpha,
            nlambda,
            lambda_min_ratio = lambda_min_ratio
        )
        call$lambda_min_ratio <- lambda_min_ratio
    }

    ## Ensure lambda is sorted in an increasing direction
    lambda <- sort(lambda)

    ## Setup cluster
    cluster <- setupCluster(
        ncores,
        cl,
        export = c("x", "y"),
        eval = {
            library(pense)
        }
    )

    ## Function returning the estimates and residuals at every lambda value
    ## (needs x and y available in the environment)
    pense_est_job <- function(segment, ...) {
        if (length(segment) == 0L) {
            return(.pense_est_pred(
                x_train = x,
                y_train = y,
                y_test = numeric(0),
                ...
            ))
        }

        return(.pense_est_pred(
            x_train = x[-segment, , drop = FALSE],
            y_train = y[-segment],
            x_test = x[segment, , drop = FALSE],
            y_test = y[segment],
            ...
        ))
    }

    ##
    ## Always do a warm-0 start to get a local optimum at each lambda value
    ##
    short_pense_options <- options
    short_pense_options$maxit <- 10L

    warm0res <- rev(pense_full(
        x = x,
        y = y,
        alpha = alpha,
        lambda_grid = rev(lambda),
        refine_it = 1L,
        nkeep = 1L,
        standardize = standardize,
        pense_options = short_pense_options,
        en_options = en_options,
        warn = FALSE
    ))

    initial_ests <- vector("list", length(warm0res))
    for (i in seq_along(initial_ests)) {
        initial_ests[[i]] <- list(
            list(
                intercept = warm0res[[i]]$intercept,
                beta = warm0res[[i]]$beta
            )
        )
    }

    ##
    ## Compute `warm_reset` "cold" initial estimators for each
    ## CV split and the full data
    ##
    get_cold_est <- function (job, ...) {
        if (length(job$segment) == 0L) {
            x_train <- x
            y_train <- y
        } else {
            x_train <- x[-job$segment, , drop = FALSE]
            y_train <- y[-job$segment]
        }

        pense_init_cold(
            x = x_train,
            y = y_train,
            lambda = job$lambda,
            ...
        )
    }

    lambda_cold_ind <- floor(seq(1, nlambda, length.out = warm_reset))
    lambda_cold <- lambda[lambda_cold_ind]
    jobs_cold_est <- lapply(lambda_cold, function (l) {
        list(
            lambda = l,
            segment = integer(0L)
        )
    })

    cv_segments <- list(integer(0L))

    if((nlambda > 1L) && (cv_k > 1L))  {
        # Create CV segments
        cv_segments <- split(
            seq_len(dx[1L]),
            sample(rep_len(seq_len(cv_k), dx[1L]))
        )

        jobs_cold_cv <- unlist(lapply(cv_segments, function (seg) {
            lapply(lambda_cold, function (l) {
                list(
                    lambda = l,
                    segment = seg
                )
            })
        }), recursive = FALSE, use.names = FALSE)

        jobs_cold_est <- c(jobs_cold_est, jobs_cold_cv)
    }

    # Compute all cold initial estimators
    tryCatch({
        cold_inits <- cluster$lapply(
            jobs_cold_est,
            get_cold_est,
            alpha = alpha,
            standardize = standardize,
            pense_options = options,
            initest_options = init_options,
            en_options = en_options
        )
    }, error = function(e) {
        cluster$stopCluster()
        stop(e)
    })

    # Group the initial estimators by lambda
    cold_inits <- split(
        cold_inits,
        rep.int(seq_len(warm_reset), length(cold_inits) %/% warm_reset)
    )

    # Prepare initial estimate list
    for (i in seq_along(lambda_cold_ind)) {
        initial_ests[[lambda_cold_ind[[i]]]] <- c(
            initial_ests[[lambda_cold_ind[[i]]]],
            unlist(cold_inits[[i]], recursive = FALSE, use.names = FALSE)
        )
    }

    ## Perform actual CV (if we have more than a single lambda value)
    if((nlambda > 1L) && (cv_k > 1L))  {
        tryCatch({
            cv_results <- cluster$lapply(
                cv_segments,
                pense_est_job,
                lambda = lambda,
                initial_ests = initial_ests,
                nkeep = init_options$keepSolutions,
                refine_it = init_options$maxitPenseRefinement,
                alpha = alpha,
                standardize = standardize,
                en_correction = !options$naiveEn,
                pense_options = options,
                en_options = en_options
            )
        }, error = function(e) {
            cluster$stopCluster()
            stop(e)
        })

        # Add CV estimates to initial estimate list
        initial_ests <- lapply(seq_along(lambda), function (j) {
            cv_ests <- lapply(cv_results, function (cvr) {
                list(
                    intercept = cvr$intercept[[j]],
                    beta = cvr$beta[ , j, drop = FALSE]
                )
            })
            c(cv_ests, initial_ests[[j]])
        })

        cv_objective_fun <- if (missing(cv_objective) || is.null(cv_objective)) {
            # robust version of the RMSPE (i.e., also taking bias into account)
            function(x) {
                st2 <- scaleTau2(x, mu.too = TRUE)
                sqrt(sum(st2^2))
            }
        } else {
            match.fun(cv_objective)
        }

        # Collect all prediction errors for each lambda sub-grid and determine
        # the optimal lambda
        all_cv_resids <- do.call(rbind, lapply(cv_results, '[[', 'residuals'))
        cv_obj <- apply(all_cv_resids, 2, function (r) {
            cv_objective_fun(r[is.finite(r)])
        })
        cv_resid_size <- apply(all_cv_resids, 2, function (r) {
            .Call(C_tau_size, r[is.finite(r)])
        })

        cv_scales <- apply(all_cv_resids, 2, function (r) {
            mscale(
                r - mean(r),
                delta = options$bdp,
                rho = "bisquare",
                cc = options$cc,
                eps = options$mscaleEps,
                maxit = options$mscaleMaxit
            )
        })

        cv_stats <- array(
            unlist(lapply(cv_results, '[[', 'sol_stats')),
            dim = c(4L, nlambda, cv_k),
            dimnames = list(
                c("obj_fun", "fold_s_scale", "beta_L1", "beta_L2"),
                NULL,
                NULL
            )
        )
        cv_stats <- apply(cv_stats, 1L, rowMeans, na.rm = TRUE)

        cv_lambda_grid <- data.frame(
            lambda = lambda,
            cvavg = cv_obj,
            resid_size = cv_resid_size,
            s_scale = cv_scales,
            cv_stats
        )
        lambda_opt <- lambda[which.min(cv_obj)]
    } else {
        cv_lambda_grid <- NULL
        lambda_opt <- lambda[1L]
    }

    cluster$stopCluster()

    ##
    ## Now compute the estimator on the full data set
    ## (using the warm-0 and the CV coefficients as additional
    ## initial estimators)
    ##

    full_results <- pense_est_job(
        integer(0L),
        lambda = lambda,
        initial_ests = initial_ests,
        nkeep = init_options$keepSolutions,
        refine_it = init_options$maxitPenseRefinement,
        alpha = alpha,
        standardize = standardize,
        en_correction = !options$naiveEn,
        pense_options = options,
        en_options = en_options
    )

    coef_ests <- rbind(
        full_results$intercept,
        full_results$beta
    )

    return(structure(list(
        residuals = full_results$residuals,
        coefficients = nameCoefVec(coef_ests, x),
        adjusted = full_results$adjusted,
        lambda = lambda,
        scale = full_results$sol_stats["scale", ],
        objective = full_results$sol_stats["objF", ],
        cv_lambda_grid = cv_lambda_grid,
        lambda_opt = lambda_opt,
        alpha = alpha,
        standardize = standardize,
        pense_options = options,
        initest_options = init_options,
        en_options = en_options,
        call = call
    ), class = "pense"))
}


##
## @param x_train,y_train training data
## @param x_test,y_test test data
## @param alpha,lambda penalty parameters
## @param initial_ests initial estimates for every lambda in the grid
## @param pense_options
## @param ... further arguments passed to pense_full
.pense_est_pred <- function(
    x_train,
    y_train,
    x_test,
    y_test,
    alpha,
    lambda,
    initial_ests,
    en_correction,
    standardize,
    ...
) {
    std_train_data <- standardize_data(x_train, y_train, standardize)
    lambda <- lambda / max(std_train_data$scale_x)
    initial_ests <- lapply(initial_ests, function(x) {
        lapply(x, std_train_data$standardize_coefs)
    })

    # Traverse from left to right (increasing lambda)
    est_all_right <- pense_full(
        x = std_train_data$xs,
        y = std_train_data$yc,
        alpha = alpha,
        lambda_grid = lambda,
        initial_ests = initial_ests,
        standardize = FALSE,
        ...
    )

    # Traverse from right to left (decreasing lambda)
    est_all_left <- rev(pense_full(
        x = std_train_data$xs,
        y = std_train_data$yc,
        alpha = alpha,
        lambda_grid = rev(lambda),
        initial_ests = rev(initial_ests),
        standardize = FALSE,
        ...
    ))

    est_all <- mapply(function (left, right) {
        if (is.null(left)) {
            return(right)
        }
        if (is.null(right)) {
            return(left)
        }

        if (isTRUE(left$objF < right$objF)) {
            return(left)
        }

        return(right)
    }, left = est_all_left, right = est_all_right, SIMPLIFY = FALSE)

    residuals <- NULL
    weights <- NULL

    ## Unstandardize coefficients
    est_all <- lapply(est_all, std_train_data$unstandardize_coefs)

    coef <- list(
        intercept = unlist(lapply(est_all, "[[", "intercept")),
        beta = do.call(cbind, lapply(est_all, "[[", "beta"))
    )

    adj_est <- NULL
    if (isTRUE(en_correction)) {
        # adj_facts needs to be on the *standardized* lambda
        adj_facts <- sqrt(1 + (1 - alpha) * lambda)

        adj_est <- mapply(function (est, adj_fact) {
            beta <- est$beta * adj_fact
            residuals <- drop(y_train - x_train %*% beta)
            intercept <- weighted.mean(residuals, est$weights)

            return(list(
                beta = beta,
                intercept = intercept,
                adj_fact = adj_fact,
                residuals = residuals - intercept
            ))
        }, est_all, adj_facts, SIMPLIFY = FALSE)
    }

    adjusted <- NULL

    if (length(y_test) == 0L) {
        residuals <- vapply(
            est_all,
            '[[',
            'residuals',
            FUN.VALUE = y_train,
            USE.NAMES = FALSE
        )

        if (isTRUE(en_correction)) {
            adjusted <- list(
                factor = unlist(lapply(adj_est, "[[", "adj_fact")),
                intercept = unlist(lapply(adj_est, "[[", "intercept"))
            )
        }
    } else {
        est_for_resid <- if (isTRUE(en_correction)) {
            adj_est
        } else {
            est_all
        }
        residuals <- vapply(
            est_for_resid,
            function(est) {
                drop(y_test - est$intercept - x_test %*% est$beta)
            },
            FUN.VALUE = y_test,
            USE.NAMES = FALSE
        )
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
        FUN.VALUE = numeric(4L),
        USE.NAMES = TRUE
    )

    return(list(
        residuals = residuals,
        adjusted = adjusted,
        intercept = coef$intercept,
        beta = coef$beta,
        sol_stats = sol_stats
    ))
}
