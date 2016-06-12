#' Perform an M-step on the PENSE result
#'
#' Performs an M-step on the S-estimator returned from \code{\link{pense}}.
#'
#' @param penseobj an object returned from a call to \code{\link{pense}}.
#' @param complete.grid should the optimal lambda be chosen from a grid of values
#'      or approximate it from the optimal lambda
#' @param version Which M-step should be performed? The newly derived ("new") or the one
#'      from mmlasso?
#' @param cv.k perform k-fold CV to choose optimal lambda (only used if
#'      \code{complete.grid = FALSE}).
#' @param nlambda the number of lambda values to try.
#' @param ncores,cl use this many cores or the supplied cluster for choosing the
#'      optima lambda. See \code{\link{pense}} for more details.
#' @importFrom robustbase .Mchi
#' @importFrom stats mad median
#' @export
mstep <- function(penseobj, complete.grid = TRUE, cv.k = 5L, nlambda = 30L,
                  ncores = getOption("mc.cores", 1L), cl = NULL) {
    X <- data.matrix(eval(penseobj$call$X, envir = parent.frame()))
    y <- eval(penseobj$call$y, envir = parent.frame())

    dX <- dim(X)

    control <- penseobj$control

    scale.x <- 1
    mux <- 0
    muy <- 0

    pense.coef <- penseobj$coefficients
    pense.lambda.opt <- penseobj$lambda.opt

    ## Standardize data and coefficients
    if (penseobj$standardize == TRUE) {
        scale.x <- apply(X, 2, mad)
        mux <- apply(X, 2, median)
        muy <- median(y)

        Xs <- scale(X, center = mux, scale = scale.x)
        yc <- y - muy

        pense.coef[1L] <- pense.coef[1L] - muy + drop(pense.coef[-1L] %*% mux)
        pense.coef[-1L] <- pense.coef[-1L] * max(scale.x)

        pense.lambda.opt <- pense.lambda.opt / max(scale.x)
    } else {
        Xs <- X
        yc <- y
    }

    ##
    ## Adjust initial scale according to Maronna & Yohai (2010)
    ## and correct for bias
    ##
    int.rho.fun <- switch (control$mscale.rho.fun,
        huber = 0L,
        bisquare = 1L,
        1L
    )
    edf <- sum(abs(pense.coef) > .Machine$double.eps) # includes intercept, no need to + 1!
    delta.adj <- control$mscale.delta * (1 - edf / dX[1L])
    delta.adj <- max(0.25, delta.adj) # Don't go below 0.25

    c.scale <- consistency.rho(delta.adj, int.rho.fun)
    scale.init <- mscale(penseobj$residuals, b = delta.adj, cc = c.scale,
                         eps = control$mscale.tol, max.it = control$mscale.maxit)

    resid.scaled <- penseobj$residuals / scale.init

    corr.fact.a <- mean(.Mchi(resid.scaled, c.scale, int.rho.fun, deriv = 1L)^2)
    corr.fact.b <- mean(.Mchi(resid.scaled, c.scale, int.rho.fun, deriv = 2L))
    corr.fact.c <- mean(.Mchi(resid.scaled, c.scale, int.rho.fun, deriv = 1L) * resid.scaled)
    scale.corr.fact <- 1 + edf / (2 * dX[1L]) * (corr.fact.a / (corr.fact.b * corr.fact.c))

    scale.init.corr <- scale.init * scale.corr.fact

    ##
    ## Select tuning constant for M-step
    ##
    c0 <- 3.44
    if (edf / dX[1L] > 0.1) {
        c0 <- 3.89
    }

    ##
    ## Adjust lambda for the M step
    ##
    lambda.opt.ls <- pense.lambda.opt * control$mscale.cc^2 / facon(control$mscale.delta)
    lambda.opt.mm <- lambda.opt.ls * 3 / c0^2


    lambda.opt.mm.cv <- lambda.opt.mm
    lambda.grid.m <- matrix(c(lambda.opt.mm, NA_real_), ncol = 2L)

    ##
    ## Optionally choose optimal lambda from a grid of candidate values
    ##
    if (isTRUE(complete.grid)) {
        ## Setup cluster
        cluster <- setupCluster(ncores, cl,
                                export = c("Xs", "yc", "version"),
                                eval = {
                                    library(pense)
                                })

        ##
        ## Find maximum lambda, starting from adjusted "optimum"
        ## in log-2 steps (i.e., always double the previous value)
        ##
        lambda.max <- lambda.opt.mm

        checkAllZero <- function(lambda, init.scale, init.coef, c0, alpha, control) {
            msres <- pensemstep(Xs, yc, cc = c0,
                                init.scale = init.scale,
                                init.coef = init.coef,
                                alpha = alpha,
                                lambda = lambda,
                                control)
            return(all(msres$beta == 0))
        }

        repeat {
            check.lambdas <- c(lambda.max, 2 * seq_len(cluster$ncores - 1L) * lambda.max)

            all.zero <- cluster$lapply(check.lambdas, checkAllZero,
                                        init.scale = scale.init.corr,
                                        init.coef = pense.coef,
                                        c0 = c0,
                                        alpha = penseobj$alpha,
                                        control = control)

            all.zero <- which(unlist(all.zero))

            if (length(all.zero) > 0) {
                lambda.max <- check.lambdas[all.zero[1L]]
                break
            }

            lambda.max <- 2 * check.lambdas[cluster$ncores]
        }

        ##
        ## Create grid of nlambda values in [lambda.min; lambda.max]
        ## where lambda.min is the adjusted minimum of the grid for the S-estimator
        ##
        lambda.tmp.grid <- lambda.grid(Xs, yc, penseobj$alpha, 2L, control, standardize = FALSE,
                                       lambda.min.ratio = eval(penseobj$call$lambda.min.ratio))

        ## Make the value even smaller than with the usual adjustment of 3 / c0^2
        ## to make sure we cover enough ground
        lambda.min <- lambda.tmp.grid[1L] / c0^2

        lambda.grid.m <- exp(seq(from = log(lambda.min), to = log(lambda.max), length = nlambda))
        lambda.grid.m <- sort(c(lambda.grid.m, lambda.opt.mm))

        ## Determine CV segments
        cv.segments <- if(cv.k > 1L) {
            split(seq_len(dX[1L]), sample(rep_len(seq_len(cv.k), dX[1L])))
        } else {
            list(integer(0))
        }

        ## Combine CV segments and lambda grid to a list of jobs
        jobgrid <- lapply(cv.segments, function(cvs) {
            lapply(lambda.grid.m, function(lvs) {
                list(segment = cvs,
                     lambda = lvs)
            })
        })

        jobgrid <- unlist(jobgrid, recursive = FALSE, use.names = FALSE)

        ## Run all jobs (combination of all CV segments and all lambda values)
        dojobcv <- function(job, init.scale, init.coef, c0, alpha, control) {
            if (length(job$segment) == 0L) {
                X.train <- Xs
                y.train <- yc
                X.test <- Xs
                y.test <- yc
            } else {
                X.train <- Xs[-job$segment, , drop = FALSE]
                y.train <- yc[-job$segment]
                X.test <- Xs[job$segment, , drop = FALSE]
                y.test <- yc[job$segment]
            }

            msres <- pense:::pensemstep(X.train, y.train, c0, init.scale, init.coef,
                                        alpha = alpha, lambda = job$lambda,
                                        control)

            return(drop(y.test - msres$intercept - X.test %*% msres$beta))
        }

        tryCatch({
            prediction.errors <- cluster$lapply(jobgrid, dojobcv,
                                                init.scale = scale.init.corr,
                                                init.coef = pense.coef,
                                                c0 = c0,
                                                alpha = penseobj$alpha,
                                                control = control)
        },
        finally = {
            cluster$stopCluster()
        })

        ## Collect all prediction errors for each lambda sub-grid and determine the optimal lambda
        cv.obj.fun <- match.fun(control$cv.objective)
        cv.performance <- unlist(lapply(split(prediction.errors, rep.int(seq_along(lambda.grid.m), cv.k)),
                                        function(preds) {
                                            preds <- unlist(preds)
                                            cv.obj.fun(preds)
                                        }), use.names = FALSE)

        lambda.opt.mm.cv <- lambda.grid.m[which.min(cv.performance)]
        lambda.grid.m <- cbind(lambda.grid.m * max(scale.x), cv.performance)
    }


    ##
    ## Compute M-estimator for lambda.opt.mm.cv
    ##
    msres <- pensemstep(Xs, yc, c0, scale.init.corr, pense.coef,
                        alpha = penseobj$alpha, lambda = lambda.opt.mm.cv,
                        control)

    ## Un-standardize the coefficients
    if (penseobj$standardize == TRUE) {
        msres$beta <- msres$beta / scale.x
        msres$intercept <- msres$intercept + muy - drop(msres$beta %*% mux)

        lambda.opt.mm.cv <- lambda.opt.mm.cv * max(scale.x)
    }

    return(structure(list(
        residuals = msres$resid,
        coefficients = nameCoefVec(c(msres$intercept, msres$beta), X),
        sest.coefficients = penseobj$coefficients,
        objective.s = penseobj$objF,
        objective = msres$objF,
        delta.adj = delta.adj,
        lambda.opt = lambda.opt.mm.cv,
        lambda.grid = lambda.grid.m,
        lambda.opt.s = penseobj$lambda.opt,
        scale = penseobj$scale,
        scale.adj = scale.init.corr,
        alpha = penseobj$alpha,
        standardize = penseobj$standardize,
        control = control,
        call = penseobj$call
    ), class = c("pensem", "pense")))
}


##
#' @useDynLib pense C_augtrans C_pen_mstep
#' @importFrom robustbase Mchi
pensemstep <- function(X, y, cc, init.scale, init.coef, alpha, lambda, control) {
    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X, dX[1L], dX[2L])
    dX[2L] <- dX[2L] + 1L

    cctrl <- initest.control(
        lambda = lambda,
        alpha = alpha,
        numIt = control$pense.maxit,
        eps = control$pense.tol^2,
        mscale.delta = control$mscale.delta,
        mscale.cc = cc,
        enpy.control = enpy.control(
            en.maxit = control$pense.en.maxit,
            en.tol = control$pense.en.tol,
            en.centering = TRUE,
            mscale.maxit = control$mscale.maxit,
            mscale.tol = control$mscale.tol,
            mscale.rho.fun = "bisquare",
            en.algorithm = control$en.algorithm
        ),

        # Not needed, set for completeness
        resid.clean.method = "proportion",
        resid.threshold = 0.5,
        resid.proportion = 0.5,
        psc.proportion = 0.5
    )

    res <- .Call(C_pen_mstep, Xtr, y, dX[1L], dX[2L], init.coef, init.scale, cctrl)

    ret <- list(
        intercept = res[[1L]][1L],
        beta = res[[1L]][-1L],
        resid = res[[2L]],
        rel.change = sqrt(res[[3L]]),
        iterations = res[[4L]],
        objF = NA_real_
    )

    ret$objF <- sum(Mchi(drop(ret$resid / init.scale), cc = cc, psi = control$mscale.rho.fun)) +
        lambda * (((1 - alpha) / 2) * sum(ret$beta^2) + alpha * sum(abs(ret$beta)))

    if (ret$rel.change > control$pense.tol) {
        warning(sprintf("M-step did not converge for lambda = %.3f", lambda))
    }

    return(ret)
}


##
#' @importFrom robustbase Mwgt MrhoInf Mchi
pensemstep.rimpl <- function(X, y, cc, init.scale, init.coef, alpha, lambda, control) {
    dX <- dim(X)
    p <- dX[2L]
    n <- dX[1L]

    current.coefs <- init.coef
    prev.coefs <- NULL

    resid <- as.vector(y - current.coefs[1L] - X %*% current.coefs[-1L])

    tol <- control$pense.tol^2
    it <- 0L
    rel.change <- Inf

    repeat {
        it <- it + 1L

        resid.scaled <- resid / init.scale
        # Mwgt is safer then Mchi in the case 0/0!
        weights <- Mwgt(resid.scaled, cc = cc, psi = control$mscale.rho.fun)

        weighted.en <- .elnet.wfit(X, y, weights = weights, alpha = alpha, lambda = lambda,
                                   centering = TRUE, maxit = control$pense.en.maxit,
                                   eps = control$pense.en.tol, warmCoef = current.coefs,
                                   en.algorithm = control$en.algorithm)

        prev.coefs <- current.coefs
        current.coefs <- weighted.en$coefficients

        resid <- weighted.en$residuals

        rel.change <- sum((prev.coefs - current.coefs)^2) / sum(prev.coefs^2)

        if (is.nan(rel.change)) { # if 0/0
            rel.change <- 0
        }

        ## Check convergence
        if (it >= control$pense.maxit || all(rel.change < tol)) {
            break
        }
    }

    ret <- list(
        intercept = current.coefs[1L],
        beta = current.coefs[-1L],
        resid = resid,
        rel.change = sqrt(rel.change),
        iterations = it,
        objF = NA_real_
    )

    ret$objF <- sum(Mchi(drop(ret$resid / init.scale), cc = cc, psi = control$mscale.rho.fun)) +
        lambda * (((1 - alpha) / 2) * sum(ret$beta^2) + alpha * sum(abs(ret$beta)))

    if (ret$rel.change > control$pense.tol) {
        warning(sprintf("M-step did not converge for lambda = %.3f", lambda))
    }

    return(ret)
}


