#' Perform an M-step after the EN S-Estimator
#'
#' Performs an M-step on the S-estimator returned from \code{\link{pense}}.
#'
#' @param penseobj an object returned from a call to \code{\link{pense}}.
#' @param lambda regularization of the initial S-estimate. Defaults to
#'      \code{penseobj$lambda_opt}, the lambda with minimal CV error.
#' @param complete_grid should the optimal lambda be chosen from a grid of
#'      values or approximate it from the optimal lambda?
#' @param cv_k perform k-fold CV to choose optimal lambda (only used if
#'      \code{complete_grid = TRUE}).
#' @param scale what initial scale to use. \code{"s"} uses the estimated
#'      M-scale of the S-regression residuals, \code{"cv"} uses the
#'      M-scale of the S-regression out-of-sample residuals. Can also
#'      be a number directly specifiying the scale to use.
#' @param nlambda the number of lambda values to try.
#' @param ncores,cl use this many cores or the supplied cluster for choosing the
#'      optima lambda. See \code{\link{pense}} for more details.
#' @param adjust_bdp should the breakdown point be adjusted based on the
#'      effective degrees of freedom?
#' @param options additional options for the M-step.
#' @param X,y,cv_objective,en_options override arguments
#'      provided to the original call to \code{\link{pense}}.
#' @importFrom robustbase .Mchi scaleTau2
#' @importFrom stats mad median
#' @importFrom Matrix drop
#' @export
mstep <- function(penseobj, lambda, complete_grid = TRUE, cv_k = 5L,
                  nlambda = 30L,
                  scale = c("s", "cv"),
                  ncores = getOption("mc.cores", 1L), cl = NULL,
                  options = mstep_options(),
                  X, y, cv_objective, en_options) {

    if (missing(X)) {
        X <- data.matrix(eval(penseobj$call$X, envir = parent.frame()))
    }

    if (missing(y)) {
        y <- eval(penseobj$call$y, envir = parent.frame())
    }

    if (missing(en_options)) {
        en_options <- penseobj$en_options
    }

    if (missing(lambda)) {
        lambda <- penseobj$lambda_opt
    }

    lambda <- .check_arg(lambda, "numeric", range = 0)
    lambda_ind <- which.min(abs(lambda - penseobj$lambda))
    if (abs(lambda - penseobj$lambda[lambda_ind]) > sqrt(.Machine$double.eps)) {
        warning("S-estimator was not computed for the given lambda. ",
                "Using the closest lambda as crude approximation.")
    }

    dX <- dim(X)

    ## store the call
    call <- match.call()
    call[[1L]] <- as.name("mstep")

    pense_lambda_opt <- penseobj$lambda_opt
    pense_coef <- coef(penseobj, lambda = lambda, exact = TRUE, sparse = TRUE)
    pense_int <- pense_coef[1L]
    pense_beta <- pense_coef[-1L, , drop = FALSE]

    residuals <- residuals(penseobj, lambda = lambda, exact = TRUE)

    ## Standardize data and coefficients
    std_data <- standardize_data(X, y, penseobj$standardize)
    Xs <- std_data$xs
    yc <- std_data$yc

    if (isTRUE(penseobj$standardize)) {
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
        lambda2 <- pense_lambda_opt * (1 - penseobj$alpha) * dX[1L]
        xtx <- crossprod(Xs[, active_set, drop = FALSE])
        hmat <- solve(xtx + lambda2 * diag(length(active_set)), xtx)
        sum(diag(hmat))
    } else {
        length(active_set)
    }

    edf <- edf + 1 # add intercept

    bdp_adj <- penseobj$pense_options$bdp
    cc_scale <- penseobj$pense_options$cc

    if (is.character(scale)) {
        scale_init <- switch(
            match.arg(scale),
            s = penseobj$scale[lambda_ind],
            cv = penseobj$cv_lambda_grid$s_scale[lambda_ind]
        )
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
        scale_corr_fact <- if (edf / dX[1L] > 0.5) {
            # This corresponds to q_T in Maronna & Yohai (2010)
            corr_fact_a <- mean(.Mchi(resid_scaled, cc_scale, 1L, deriv = 1L)^2)
            corr_fact_b <- mean(.Mchi(resid_scaled, cc_scale, 1L, deriv = 2L))
            corr_fact_c <- mean(.Mchi(resid_scaled, cc_scale, 1L, deriv = 1L) *
                                    resid_scaled)

            1 + edf / (2 * dX[1L]) * (corr_fact_a / (corr_fact_b * corr_fact_c))
        } else if (edf / dX[1L] > 0.1) {
            ## this is q_E in Maronna & Yohai (2010)
            1 / (1 - (1.29 - 6.02 / dX[1L]) * edf / dX[1L])
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
    options$cc <- if (edf / dX[1L] > 0.33) {
        4.2
    } else if (edf / dX[1L] > 0.2) {
        4
    } else if (edf / dX[1L] > 0.1) {
        3.7
    } else {
        3.44
    }

    ##
    ## Adjust lambda for the M step
    ##
    lambda_opt_m <- pense_lambda_opt * 1.5 / options$cc^2

    lambda_opt_m_cv <- lambda_opt_m
    lambda_grid_m <- data.frame(
        lambda = lambda_opt_m_cv,
        cvavg = NA_real_,
        cvsd = NA_real_
    )

    ##
    ## Optionally choose optimal lambda from a grid of candidate values
    ##
    if (isTRUE(complete_grid)) {
        ## Setup cluster
        cluster <- setupCluster(
            ncores,
            cl,
            export = c("Xs", "yc"),
            eval = {
                library(pense)
            }
        )

        ##
        ## Find maximum lambda, starting from adjusted "optimum"
        ## in log-2 steps (i.e., always double the previous value)
        ##
        lambda_max <- lambda_min <- lambda_opt_m

        get_coef_norm <- function(...) {
            msres <- pensemstep(
                Xs,
                yc,
                ...
            )
            # msres$beta is of type `dgCMatrix`
            return(sum(abs(msres$beta@x)))
        }

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
                alpha = penseobj$alpha,
                options = options,
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
                    alpha = penseobj$alpha,
                    options = options,
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

        # We have found a good lambda_max, now let's look for a good lambda_min
        lambda_min_break <- eval(penseobj$call$lambda_min_ratio) * lambda_max
        lambda_min_break <- ifelse(is.null(lambda_min_break), 1e-5 * lambda_max,
                                   lambda_min_break)

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
                    alpha = penseobj$alpha,
                    options = options,
                    en_options = en_options
                )

                coef_norm <- unlist(coef_norm)
                coef_norm_diff <- coef_norm[-length(coef_norm)] / coef_norm[-1L]
                const_norm <- which(abs(coef_norm_diff - 1) < options$eps)
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

        ## Determine CV segments
        cv_segments <- if(cv_k > 1L) {
            split(seq_len(dX[1L]), sample(rep_len(seq_len(cv_k), dX[1L])))
        } else {
            list(integer(0))
        }

        ## Combine CV segments and lambda grid to a list of jobs
        jobgrid <- lapply(cv_segments, function(cvs) {
            lapply(lambda_grid_m, function(lvs) {
                list(
                    segment = cvs,
                    lambda = lvs
                )
            })
        })

        jobgrid <- unlist(jobgrid, recursive = FALSE, use.names = FALSE)

        ## Run all jobs (combination of all CV segments and all lambda values)
        dojobcv <- function(job, ...) {
            if (length(job$segment) == 0L) {
                X_train <- Xs
                y_train <- yc
                X_test <- Xs
                y_test <- yc
            } else {
                X_train <- Xs[-job$segment, , drop = FALSE]
                y_train <- yc[-job$segment]
                X_test <- Xs[job$segment, , drop = FALSE]
                y_test <- yc[job$segment]
            }

            msres <- pensemstep(
                X_train,
                y_train,
                lambda = job$lambda,
                ...
            )

            return(drop(y_test - msres$intercept - X_test %*% msres$beta))
        }

        tryCatch({
            pred_errors <- cluster$lapply(
                jobgrid,
                dojobcv,
                alpha = penseobj$alpha,
                init_scale = scale_init_corr,
                init_int = pense_int,
                init_coef = pense_beta,
                options = options,
                en_options = en_options
            )
        },
        finally = {
            cluster$stopCluster()
        })

        ## Collect all prediction errors for each lambda sub-grid and
        ## determine the optimal lambda
        cv_objective_fun <- if (missing (cv_objective)) {
            if (is.null(penseobj$call$cv_objective)) {
                scaleTau2
            } else {
                match.fun(penseobj$call$cv_objective)
            }
        } else {
                match.fun(cv_objective)
        }

        cv_obj <- unlist(
            lapply(
                split(pred_errors, rep.int(seq_along(lambda_grid_m), cv_k)),
                function(preds) {
                    pred_resids <- unlist(preds)
                    cv_objective_fun(pred_resids[is.finite(pred_resids)])
                }
            )
        )

        lambda_opt_m_cv <- lambda_grid_m[which.min(cv_obj)]
        lambda_grid_m <- data.frame(
            lambda = lambda_grid_m * max(std_data$scale_x),
            cvavg = cv_obj
        )
    }

    ##
    ## Compute M-estimator for lambda.opt.mm.cv
    ##
    msres <- pensemstep(
        Xs,
        yc,
        init_scale = scale_init_corr,
        init_int = pense_int,
        init_coef = pense_beta,
        alpha = penseobj$alpha,
        lambda = lambda_opt_m_cv,
        options = options,
        en_options = en_options
    )

    ## Un-standardize the coefficients
    if (isTRUE(penseobj$standardize)) {
        msres <- std_data$unstandardize_coefs(msres)
        lambda_opt_m_cv <- lambda_opt_m_cv * max(std_data$scale_x)
    }

    return(structure(list(
        residuals = msres$resid,
        coefficients = nameCoefVec(rbind(msres$intercept, msres$beta), X),
        objective = msres$objF,
        bdp = bdp_adj,
        lambda_opt = lambda_opt_m_cv,
        lambda = lambda_opt_m_cv,
        cv_lambda_grid = lambda_grid_m,
        scale = scale_init_corr,
        alpha = penseobj$alpha,
        standardize = penseobj$standardize,
        sest = penseobj,
        call = call
    ), class = c("pensem", "pense")))
}

##
#' @useDynLib pense C_augtrans C_pen_mstep_sp
#' @importFrom robustbase Mchi
#' @importFrom methods is
#' @importFrom Matrix Matrix
#' @importClassesFrom Matrix dgCMatrix
pensemstep <- function(X, y, init_scale, init_int, init_coef, alpha, lambda,
                       options, en_options) {
    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X)
    dX[2L] <- dX[2L] + 1L

    if (!is(init_coef, "dgCMatrix")) {
        init_coef <- Matrix(init_coef, sparse = TRUE, ncol = 1L)
    }

    res <- .Call(C_pen_mstep_sp, Xtr, y, init_int, init_coef, init_scale,
                 alpha, lambda, options, en_options)

    res$objF <- sum(Mchi(drop(res$residuals / init_scale),
                         cc = options$cc, psi = "bisquare")) +
        lambda * (
            0.5 * (1 - alpha) * sum(res$beta^2) + alpha * sum(abs(res$beta))
        )

    ##
    ## Check if the M-step converged.
    ## Be extra careful with the comparison as rel.change may be NaN or NA
    ##
    if (!isTRUE(res$rel_change < options$eps)) {
        warning(sprintf("M-step did not converge for lambda = %.6f", lambda))
    }

    return(res)
}
