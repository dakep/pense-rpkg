#' Perform an M-step after the EN S-Estimator
#'
#' Compute the PENSEM estimate, an efficient and robust
#' elastic net estimator for linear regression.
#'
#' Performs an M-step using the S-estimator at the optimal penalty
#' parameter as returned from \code{\link{pense}} as the initial
#' estimate. For "fat" datasets, the initial scale as returned by the
#' S-estimate is adjusted according to Maronna & Yohai (2010).
#'
#' @param x either a numeric data matrix or a fitted PENSE estimate obtained
#'      from \code{\link{pense}}.
#' @param alpha elastic net mixing parameter with \eqn{0 \le \alpha \le 1}.
#'      \code{alpha = 1} is the LASSO penalty, and \code{alpha = 0} the Ridge
#'      penalty. If a \code{pense} object is supplied as first argument,
#' @param nlambda if \code{lambda} is not given,
#'      a grid of \code{nlambda} lambda values is generated based on the data.
#' @param lambda a single value or a grid of values for the regularization
#'      parameter of the M-step.
#'      Assumed to be on the same scale as the data.
#'      If missing, a grid of lambda
#'      values is automatically generated (see parameter \code{nlambda}).
#'      If supplied and \code{standardize = TRUE}, the lambda values will be
#'      adjusted accordingly.
#' @param lambda_min_ratio If the grid should be chosen automatically, the
#'      ratio of the smallest lambda to the (computed) largest lambda.
#' @param standardize should the data be standardized robustly? Estimates
#'      are returned on the original scale. Defaults to \code{TRUE}.
#' @param cv_k perform k-fold CV to choose the optimal lambda for prediction.
#' @param cv_objective a function (name) to compute the CV performance.
#'      By default, the robust tau-scale is used.
#' @param ncores,cl use multiple cores or the supplied cluster for the
#'      cross-validation. See \code{\link{pense}} for more details.
#' @param mm_options additional options for the M-step.
#' @param en_options additional options for the EN algorithm.
#'      See \code{\link{elnet}} and \code{\link{en_options}} for details.
#' @param ... currently ignored.
#'
#' @return An object of class \code{"pensem"}. All elements as an object
#'      of class \code{\link{pense}} as well as the following:
#' \item{init_scale}{the initial scale estimate used in the M step.}
#' \item{sest}{the PENSE estimate used to initialize the M step.}
#' \item{bdp}{breakdown point of the MM-estimator.}
#'
#' @references Maronna, R. and Yohai, V. (2010).
#'      Correcting MM estimates for "fat" data sets.
#'      \emph{Computational Statistics & Data Analysis},
#'      \bold{54}:31683173.
#'
#' @seealso \code{\link{pense}} to compute only the S-estimator.
#'
#' @importFrom robustbase .Mchi scaleTau2
#' @importFrom stats mad median coef weighted.mean
#' @importFrom Matrix drop
#'
#' @example examples/pensem-1.R
#'
#' @rdname pensem
#' @export
pensem <- function(x, ...) {
    UseMethod("pensem", x)
}

#' @param y numeric response vector.
#' @param lambda_s regularization parameter for the \emph{S-estimator}.
#'      If missing, a grid of lambda values is chosen automatically.
#'      If \code{standardize = TRUE}, the lambda values will be
#'      adjusted accordingly.
#' @param initial how to initialize the estimator at a new lambda in the grid.
#'      The default, \code{"warm"}, computes a cold initial estimator at several
#'      lambda values and uses the PENSE coefficient to warm-start the
#'      estimator at the next larger lambda value. At the largest value in
#'      the lambda grid, PENSE will be initialized with the 0-vector.
#'      \code{"cold"} computes the full initial estimator at
#'      every lambda value.
#' @param warm_reset if \code{initial = "warm"} (default), how many cold initial
#'      estimates be computed?
#' @param s_options additional options for the PENSE algorithm.
#'      See \code{\link{pense_options}} for details.
#' @param init_options additional options for computing the cold initial
#'      estimates.
#'      Ignored if \code{initial = "warm"} and \code{warm_reset = 0}.
#'      See \code{\link{initest_options}} for details.
#'
#' @rdname pensem
#' @method pensem default
#' @export
pensem.default <- function(
    x, y,
    alpha = 0.5,
    nlambda = 50,
    lambda,
    lambda_s,
    lambda_min_ratio,
    standardize = TRUE,
    initial = c("warm", "cold"),
    warm_reset = 10,
    cv_k = 5, cv_objective,
    ncores = getOption("mc.cores", 1L), cl = NULL,
    s_options = pense_options(),
    mm_options = mstep_options(),
    init_options = initest_options(),
    en_options = en_options_aug_lars(),
    ...
) {
    lambda_s <- if (!missing(lambda_s)) {
        lambda_s
    } else {
        NULL
    }

    lambda_m <- if (!missing(lambda)) {
        lambda
    } else {
        NULL
    }

    lambda_min_ratio <- if (!missing(lambda_min_ratio)) {
        lambda_min_ratio
    } else {
        NULL
    }

    cv_objective <- if (!missing(cv_objective)) {
        cv_objective
    } else {
        NULL
    }

    pense_est <- pense(
        x = x,
        y = y,
        alpha = alpha,
        nlambda = nlambda,
        lambda = lambda_s,
        lambda_min_ratio = lambda_min_ratio,
        standardize = standardize,
        initial = initial,
        warm_reset = warm_reset,
        cv_k = cv_k,
        cv_objective = cv_objective,
        ncores = ncores,
        cl = cl,
        options = s_options,
        init_options = init_options,
        en_options = en_options
    )

    return(pensem.pense(
        pense_est,
        alpha = alpha,
        nlambda = nlambda,
        lambda = lambda_m,
        lambda_min_ratio = lambda_min_ratio,
        standardize = standardize,
        cv_k = cv_k,
        cv_objective = cv_objective,
        ncores = ncores,
        cl = cl,
        mm_options = mm_options,
        en_options = en_options,
        x_train = x,
        y_train = y
    ))
}



#' Refine an already computed PENSE with an additional M-step
#'
#' @param scale initial scale estimate for the M step. By default the
#'      S-scale from the initial estimator (\code{x}) is used.
#' @param x_train,y_train override arguments
#'      provided to the original call to \code{\link{pense}}.
#'
#' @rdname pensem
#' @method pensem pense
#' @export
pensem.pense <- function(
    x,
    alpha,
    scale,
    nlambda = 50,
    lambda,
    lambda_min_ratio,
    standardize,
    cv_k = 5, cv_objective,
    ncores = getOption("mc.cores", 1L), cl = NULL,
    mm_options = mstep_options(),
    en_options,
    x_train, y_train,
    ...
) {
    ## store the call
    call <- match.call()
    call[[1L]] <- as.name("pensem")

    penseobj <- x

    if (missing(x_train) || is.null(x_train)) {
        call$x_train <- penseobj$call$x
        x <- data.matrix(eval(call$x_train, envir = parent.frame()))
    } else {
        x <- x_train
    }

    if (missing(y_train) || is.null(y_train)) {
        call$y_train <- penseobj$call$y
        y <- eval(call$y_train, envir = parent.frame())
    } else {
        y <- y_train
    }

    if (missing(alpha) || is.null(alpha)) {
        alpha <- penseobj$alpha
    }

    if (missing(en_options) || is.null(en_options)) {
        en_options <- penseobj$en_options
    }

    if (missing(standardize) || is.null(standardize)) {
        standardize <- penseobj$standardize
    }

    if (missing(lambda_min_ratio) || is.null(lambda_min_ratio)) {
        lambda_min_ratio <- if (is.null(penseobj$call$lambda_min_ratio)) {
            eval(penseobj$call$lambda_min_ratio)
        } else {
            .default_lambda_min_ratio(x)
        }
    }

    pense_lambda_opt <- penseobj$lambda_opt
    lambda_ind <- which.min(abs(pense_lambda_opt - penseobj$lambda))

    if (!missing(lambda) && !is.null(lambda)) {
        lambda <- .check_arg(lambda, "numeric", range = 0, length = NULL)
    }

    dx <- dim(x)

    pense_coef <- coef(penseobj, lambda = pense_lambda_opt, exact = TRUE,
                       sparse = TRUE, correction = FALSE)
    pense_int <- pense_coef[1L]
    pense_beta <- pense_coef[-1L, , drop = FALSE]

    residuals <- drop(y - x %*% pense_beta - pense_int)

    ## Standardize data and coefficients
    std_data <- standardize_data(x, y, standardize)
    xs <- std_data$xs
    yc <- std_data$yc

    if (isTRUE(standardize)) {
        pense_coefs <- std_data$standardize_coefs(list(
            intercept = pense_int,
            beta = pense_beta
        ))
        pense_int <- pense_coefs$intercept
        pense_beta <- pense_coefs$beta
    }

    pense_lambda_opt <- pense_lambda_opt / max(std_data$scale_x)

    ##
    ## Compute effective degrees of freedom
    ## according to
    ## Tibshirani, Ryan J.; Taylor, Jonathan. Degrees of freedom in lasso problems.
    ## Ann. Statist. 40 (2012), no. 2, 1198--1232. doi:10.1214/12-AOS1003.
    ##
    active_set <- pense_beta@i + 1L

    edf <- if ((penseobj$alpha < 1) && (length(active_set) > 0L)) {
        # this is not the lambda_2 in the objective used for optimization
        # since the optimization uses differently scaled objective
        lambda2 <- pense_lambda_opt * (1 - penseobj$alpha) * dx[1L]
        xtx <- crossprod(xs[, active_set, drop = FALSE])
        hmat <- solve(xtx + lambda2 * diag(length(active_set)), xtx)
        sum(diag(hmat))
    } else {
        length(active_set)
    }

    edf <- edf + 1 # add intercept

    bdp_adj <- penseobj$pense_options$bdp
    cc_scale <- penseobj$pense_options$cc

    if (missing(scale)) {
        scale_init <- penseobj$scale[lambda_ind]
        adjust_scale <- TRUE
    } else {
        scale_init <- .check_arg(scale, "numeric", range = 0)
        adjust_scale <- FALSE
    }

    resid_scaled <- residuals / scale_init

    ##
    ## Adjust the scale for "fat" datasets
    ##
    scale_init_corr <- if (adjust_scale) {
        scale_corr_fact <- if (edf / dx[1L] > 0.5) {
            # This corresponds to q_T in Maronna & Yohai (2010)
            corr_fact_a <- mean(.Mchi(resid_scaled, cc_scale, 1L, deriv = 1L)^2)
            corr_fact_b <- mean(.Mchi(resid_scaled, cc_scale, 1L, deriv = 2L))
            corr_fact_c <- mean(.Mchi(resid_scaled, cc_scale, 1L, deriv = 1L) *
                                    resid_scaled)

            1 + edf / (2 * dx[1L]) * (corr_fact_a / (corr_fact_b * corr_fact_c))
        } else if (edf / dx[1L] > 0.1) {
            ## this is q_E in Maronna & Yohai (2010)
            1 / (1 - (1.29 - 6.02 / dx[1L]) * edf / dx[1L])
        } else {
            1
        }

        scale_init * scale_corr_fact
    } else {
        scale_init
    }

    ##
    ## Select tuning constant for the M-step for "fat" datasets
    ##
    mm_options$cc <- if (edf / dx[1L] > 0.33) {
        4.2
    } else if (edf / dx[1L] > 0.2) {
        4
    } else if (edf / dx[1L] > 0.1) {
        3.7
    } else {
        3.44
    }

    ## Setup cluster
    cluster <- setupCluster(
        ncores,
        cl,
        export = c("xs", "yc"),
        eval = {
            library(pense)
        }
    )

    ##
    ## If no lambda was provided, choose the grid automatically
    ##
    if (missing(lambda) || is.null(lambda)) {
        ##
        ## Adjust lambda for the M step
        ##
        lambda_opt_m <- pense_lambda_opt * 1.5 / mm_options$cc^2
        lambda_opt_m_cv <- lambda_opt_m

        get_coef_norm <- function(...) {
            msres <- pensemstep(
                xs,
                yc,
                ...
            )
            # msres$beta is of type `dgCMatrix`
            return(sum(abs(msres$beta@x)))
        }

        if (alpha < 0.1) {
            lambda_max <- max(penseobj$lambda) * 1.5 /
                (mm_options$cc^2 * max(std_data$scale_x))
        } else {
            ##
            ## Find maximum lambda, starting from adjusted "optimum"
            ## in log-2 steps (i.e., always double the previous value)
            ##
            lambda_max <- lambda_min <- lambda_opt_m

            ##
            ## Get largest lambda necessary
            ##
            moved_away <- FALSE
            repeat {
                check_lambdas <- c(lambda_max, (2^seq_len(cluster$ncores - 1L)) *
                                       lambda_max)

                zero_norm <- cluster$lapply(
                    check_lambdas,
                    get_coef_norm,
                    init_scale = scale_init_corr,
                    init_int = pense_int,
                    init_coef = pense_beta,
                    alpha = alpha,
                    options = mm_options,
                    en_options = en_options
                )

                zero_norm <- which(unlist(zero_norm) < .Machine$double.eps)

                if (length(zero_norm) > 0) {
                    # Check if we moved away from the initial guess
                    moved_away <- moved_away | !identical(zero_norm[1L], 1L)
                    lambda_max <- check_lambdas[zero_norm[1L]]
                    break
                }

                lambda_max <- 2 * check_lambdas[cluster$ncores]
                moved_away <- TRUE
            }

            if (!moved_away) {
                # If we did not move away from the initial guess, we have to check
                # the next smallest step if this may lead to all-zeros as well
                # check lambda_max / 2, lambda_max / 4, lambda_max / (2^...)
                repeat {
                    check_lambdas <- lambda_max * 0.5^seq_len(cluster$ncores)

                    zero_norm <- cluster$lapply(
                        check_lambdas,
                        get_coef_norm,
                        init_scale = scale_init_corr,
                        init_int = pense_int,
                        init_coef = pense_beta,
                        alpha = alpha,
                        options = mm_options,
                        en_options = en_options
                    )

                    not_zero_norm <- which(unlist(zero_norm) > .Machine$double.eps)

                    if (length(not_zero_norm) > 0L) {
                        choose <- not_zero_norm[1L] - 1L
                        lambda_max <- ifelse(choose < 1L, lambda_max,
                                             check_lambdas[choose])
                        break
                    }

                    lambda_max <- check_lambdas[cluster$ncores] * 0.5

                    if (isTRUE(lambda_max < .Machine$double.eps)) {
                        lambda_max <- lambda_opt_m
                        warning("M-step results in the zero-vector for ",
                                "all penalty values.")
                        break
                    }
                }
            }
        }

        # We have found a good lambda_max, now let's look for a good lambda_min
        lambda_min_break <- lambda_min_ratio * lambda_max
        lambda_min_break <- ifelse(
            isTRUE(is.finite(lambda_min_break)),
            lambda_min_break,
            1e-5 * lambda_max
        )

        lambda_min <- min(lambda_max, lambda_opt_m)
        if (lambda_max > .Machine$double.eps) {
            repeat {
                check_lambdas <- lambda_min * 0.5^seq_len(max(5L, cluster$ncores))

                coef_norm <- cluster$lapply(
                    check_lambdas,
                    get_coef_norm,
                    init_scale = scale_init_corr,
                    init_int = pense_int,
                    init_coef = pense_beta,
                    alpha = alpha,
                    options = mm_options,
                    en_options = en_options
                )

                coef_norm <- unlist(coef_norm)
                coef_norm_diff <- coef_norm[-length(coef_norm)] / coef_norm[-1L]
                const_norm <- which(abs(coef_norm_diff - 1) < mm_options$eps)
                if (length(const_norm) > 0L) {
                    lambda_min <- check_lambdas[[const_norm[[1L]]]]
                    break
                }

                lambda_min <- min(check_lambdas)
                if (lambda_min < lambda_min_break) {
                    lambda_min <- lambda_min_break
                    break
                }
            }
        }

        ##
        ## Create grid of nlambda values in [lambda.min; lambda.max]
        ## where lambda.min is the adjusted minimum of the grid for the S-estimator
        ##
        lambda_grid_m <- exp(seq(
            from = log(lambda_min),
            to = log(lambda_max),
            length = nlambda
        ))
        lambda_grid_m <- sort(c(lambda_grid_m, lambda_opt_m))
    } else {
        lambda_grid_m <- sort(lambda / max(std_data$scale_x))
    }

    ##
    ## Choose optimal lambda from a grid of candidate values
    ##
    if (cv_k > 1L) {
        cv_segments <- split(
            seq_len(dx[1L]),
            sample(rep_len(seq_len(cv_k), dx[1L]))
        )

        ## Check if there are more cores than CV segments available
        if (isTRUE(cluster$ncores / cv_k > 1)) {
            lg_length <- rep.int(cluster$ncores %/% cv_k, cv_k)
            lg_length[seq_len(cluster$ncores %% cv_k)] <- lg_length[1L] + 1L
            lg_length <- pmin(
                lg_length,
                length(lambda_grid_m)
            )
        } else {
            lg_length <- rep.int(1L, max(1L, cv_k))
        }

        ## Combine CV segments and lambda grid to a list of jobs
        jobgrid <- unlist(
            mapply(function (cvs, num_lgs) {
                lambda_subgrids <- split(
                    lambda_grid_m,
                    sort(rep(seq_len(num_lgs), length = length(lambda_grid_m)))
                )

                lapply(lambda_subgrids, function(lg) {
                    list(
                        segment = cvs,
                        lambda = lg
                    )
                })
            }, cvs = cv_segments, num_lgs = lg_length, SIMPLIFY = FALSE),
            recursive = FALSE, use.names = FALSE
        )

        ## Run all jobs (combination of all CV segments and all lambda values)
        dojobcv <- function(job, ...) {
            x_train <- xs[-job$segment, , drop = FALSE]
            y_train <- yc[-job$segment]
            x_test <- xs[job$segment, , drop = FALSE]
            y_test <- yc[job$segment]

            msres <- pensemstep(
                x_train,
                y_train,
                lambda = job$lambda,
                ...
            )

            return(y_test - sweep(
                x_test %*% msres$beta,
                2,
                msres$intercept,
                check.margin = FALSE
            ))
        }

        tryCatch({
            pred_errors <- cluster$lapply(
                jobgrid,
                dojobcv,
                alpha = alpha,
                init_scale = scale_init_corr,
                init_int = pense_int,
                init_coef = pense_beta,
                options = mm_options,
                en_options = en_options
            )
        },
        finally = {
            cluster$stopCluster()
        })

        ## Collect all prediction errors for each lambda sub-grid and
        ## determine the optimal lambda
        cv_objective_fun <- if (missing(cv_objective) || is.null(cv_objective)) {
            if (is.null(penseobj$call$cv_objective) ||
                is.null(eval(penseobj$call$cv_objective))) {
                function(x) {
                    st2 <- scaleTau2(x, mu.too = TRUE)
                    sqrt(sum(st2^2))
                }
            } else {
                match.fun(penseobj$call$cv_objective)
            }
        } else {
            match.fun(cv_objective)
        }

        cv_pred_errors <- do.call(rbind, lapply(
            split(pred_errors, rep.int(seq_along(lg_length), times = lg_length)),
            function(preds) {
                do.call(cbind, preds)
            }
        ))

        cv_obj <- apply(cv_pred_errors, 2, function (r) {
            cv_objective_fun(r[is.finite(r)])
        })
        cv_resid_size <- apply(cv_pred_errors, 2, function (r) {
            .Call(C_tau_size, r[is.finite(r)])
        })

        lambda_opt_m_cv <- lambda_grid_m[which.min(cv_obj)]
        cv_lambda_grid <- data.frame(
            lambda = lambda_grid_m,
            cvavg = cv_obj,
            resid_size = cv_resid_size
        )
    } else {
        lambda_opt_m_cv <- min(lambda_grid_m)
        cv_lambda_grid <- data.frame(
            lambda = lambda_grid_m,
            cvavg = NA_real_,
            resid_size = NA_real_
        )
    }

    ##
    ## Compute M-estimator for all lambda in the grid
    ##
    msres <- pensemstep(
        xs,
        yc,
        init_scale = scale_init_corr,
        init_int = pense_int,
        init_coef = pense_beta,
        alpha = alpha,
        lambda = lambda_grid_m,
        options = mm_options,
        en_options = en_options
    )

    ## Adjustment factors are based on the *standardized* coefficients
    adj_fact <- sqrt(1 + (1 - penseobj$alpha) * lambda_grid_m)

    ## Un-standardize the coefficients
    if (isTRUE(standardize)) {
        msres <- std_data$unstandardize_coefs(msres)
        lambda_opt_m_cv <- lambda_opt_m_cv * max(std_data$scale_x)
        lambda_grid_m <- lambda_grid_m * max(std_data$scale_x)
        cv_lambda_grid$lambda <- cv_lambda_grid$lambda * max(std_data$scale_x)
    }

    ## Adjusted intercepts are computed directly for the *unstandardized* coefficients
    adjusted <- list(
        factor = adj_fact,
        intercept = unlist(lapply(seq_along(lambda_grid_m), function (i) {
            weighted.mean(
                drop(y - x %*% msres$beta[ , i, drop = FALSE] * adj_fact[i]),
                msres$weights[ , i]
            )
        }), use.names = FALSE)
    )

    return(structure(list(
        residuals = msres$residuals,
        coefficients = nameCoefVec(rbind(msres$intercept, msres$beta), x),
        adjusted = adjusted,
        objective = msres$objF,
        bdp = bdp_adj,
        lambda_opt = lambda_opt_m_cv,
        lambda = lambda_grid_m,
        cv_lambda_grid = cv_lambda_grid,
        init_scale = scale_init_corr,
        alpha = alpha,
        standardize = standardize,
        sest = penseobj,
        options = mm_options,
        en_options = en_options,
        call = call
    ), class = c("pensem", "pense")))
}


## This performs the actual M-step on the given initial estimate.
#' @useDynLib pense, .registration = TRUE
#' @importFrom robustbase Mchi
#' @importFrom methods is
#' @importFrom Matrix Matrix
#' @importClassesFrom Matrix dgCMatrix
pensemstep <- function(x, y, init_scale, init_int, init_coef, alpha, lambda,
                       options, en_options) {
    dx <- dim(x)

    xtr <- .Call(C_augtrans, x)
    dx[2L] <- dx[2L] + 1L

    if (!is(init_coef, "dgCMatrix")) {
        init_coef <- Matrix(init_coef, sparse = TRUE, ncol = 1L)
    }

    res <- .Call(C_pen_mstep_sp, xtr, y, init_int, init_coef, init_scale,
                 alpha, lambda, options, en_options)

    ##
    ## Check if the M-step converged.
    ## Be extra careful with the comparison as rel.change may be NaN or NA
    ##
    if (!isTRUE(any(res$rel_change < options$eps))) {
        not_converged <- which(res$rel_change < options$eps)
        warning(sprintf(
            "M-step did not converge for lambda = %s",
            paste(sprintf("%g", lambda[not_converged]), collapse = ", ")
        ))
    }

    return(res)
}

