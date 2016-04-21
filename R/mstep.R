#' Perform an M-step on the PENSE result
#'
#' Performs an M-step on the L1 S-estimator returned from \code{\link{pense}}.
#'
#' @param penseobj an object returned from a call to \code{\link{pense}}.
#' @param complete.grid should the optimal lambda be chosen from a grid of values
#'      or approximate it from the optimal lambda
#' @param cv.k perform k-fold CV to choose optimal lambda (only used if
#'      \code{complete.grid = FALSE}).
#' @param nlambda the number of lambda values to try.
#' @param ncores,cl use this many cores or the supplied cluster for choosing the
#'      optima lambda. See \code{\link{pense}} for more details.
#' @importFrom robustbase .Mchi
#' @importFrom stats mad median
#' @export
mstep <- function(penseobj, complete.grid, cv.k = 5L, nlambda = 30L,
                  ncores = getOption("mc.cores", 2L), cl = NULL) {
    X <- data.matrix(eval(penseobj$call$X))
    y <- eval(penseobj$call$y)

    if (penseobj$alpha != 1) {
        stop("M-step not yet implemented for alpha != 1.")
    }

    dX <- dim(X)

    control <- penseobj$control

    scale.x <- 1
    mux <- 0
    muy <- 0

    ## Standardize data and coefficients
    if (penseobj$standardize == TRUE) {
        scale.x <- apply(X, 2, mad)
        mux <- apply(X, 2, median)
        muy <- median(y)

        Xs <- scale(X, center = mux, scale = scale.x)
        yc <- y - muy

        pense.coef <- penseobj$coefficients
        pense.coef[1L] <- pense.coef[1L] - muy + drop(pense.coef[-1L] %*% mux)
        pense.coef[-1L] <- pense.coef[-1L] * max(scale.x)
    } else {
        Xs <- X
        yc <- y
        pense.coef <- penseobj$coefficients
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
    delta.adj <- control$mscale.delta * (1 - edf / n)
    delta.adj <- max(0.25, delta.adj) # Don't go below 0.25

    c.scale <- consistency.rho(delta.adj, int.rho.fun)
    scale.init <- mscale(penseobj$residuals, b = delta.adj, cc = c.scale,
                         eps = control$mscale.tol, max.it = control$mscale.maxit)

    resid.scaled <- penseobj$residuals / scale.init

    corr.fact.a <- mean(.Mchi(resid.scaled, c.scale, int.rho.fun, deriv = 1L)^2)
    corr.fact.b <- mean(.Mchi(resid.scaled, c.scale, int.rho.fun, deriv = 2L))
    corr.fact.c <- mean(.Mchi(resid.scaled, c.scale, int.rho.fun, deriv = 1L) * resid.scaled)
    scale.corr.fact <- 1 + edf / (2 * n) * (corr.fact.a / (corr.fact.b * corr.fact.c))

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
    lambda.opt.ls <- penseobj$lambda.opt * control$mscale.cc^2 / facon(control$mscale.delta)
    lambda.opt.mm <- lambda.opt.ls * 3 / c0^2 * scale.init.corr

    lambda.opt.mm.cv <- lambda.opt.mm
    lambda.grid.m <- matrix(c(lambda.opt.mm, NA_real_), ncol = 2L)

    ##
    ## Optionally choose optimal lambda from a grid of candidate values
    ##
    if (isTRUE(complete.grid)) {
        lambda.grid.m <- lambda.grid(Xs, yc, 1, nlambda, control, standardize = FALSE,
                                     lambda.min.ratio = eval(penseobj$call$lambda.min.ratio))
        lambda.grid.m <- sort(c(lambda.grid.m, lambda.opt.mm))

        ## Setup cluster
        cluster <- setupCluster(ncores, cl,
                                export = c("Xs", "yc"),
                                eval = {
                                    library(penseinit)
                                })

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
        dojobcv <- function(job, init.scale, init.coef, c0, control) {
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

            coefsm <- penseinit:::pensemstep(X.train, y.train, c0, init.scale, init.coef, job$lambda,
                                             control)

            return(drop(y.test - coefsm[1L] - X.test %*% coefsm[-1L]))
        }

        tryCatch({
            prediction.errors <- cluster$lapply(jobgrid, dojobcv,
                                                init.scale = scale.init.corr,
                                                init.coef = pense.coef,
                                                c0 = c0,
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
    coefs.mm <- pensemstep(Xs, yc, c0, scale.init.corr, pense.coef, lambda.opt.mm.cv,
                           control)

    coefs.mm <- drop(coefs.mm)

    residuals <- yc - coefs.mm[1L] - Xs %*% coefs.mm[-1L]

    ## Un-standardize the coefficients
    if (penseobj$standardize == TRUE) {
        coefs.mm[-1L] <- coefs.mm[-1L] / scale.x
        coefs.mm[1L] <- coefs.mm[1L] + muy - drop(coefs.mm[-1L] %*% mux)

        lambda.opt.mm.cv <- lambda.opt.mm.cv * max(scale.x)
    }


    return(structure(list(
        residuals = residuals,
        coefficients = nameCoefVec(coefs.mm, X),
        sest.coefficients = penseobj$coefficients,
        objective = penseobj$objF,
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
    ), class = c("pensem1", "pense")))
}


## Get the constant needed for consistency for the given delta
## and the given rho function
#' @importFrom robustbase .Mchi
consistency.rho <- function(delta, int.rho.fun) {
    integrand <- function(x, cc) {
        dnorm(x) * .Mchi(x, cc, int.rho.fun)
    }

    expectation <- function(cc, delta) {
        integrate(integrand, lower = -Inf, upper = Inf, cc)$value - delta
    }

    uniroot(expectation, interval = c(0.3, 10), delta)$root
}
