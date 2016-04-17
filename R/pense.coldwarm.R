## PENSE with a cold start at the beginning and subsequently uses the previous
## parameter estimate as warm-start
##
## @param lambda.grid A grid of lambda values NOT ADJUSTED for sample size!
#' @importFrom stats mad median
pense.coldwarm <- function(X, y, alpha, lambda.grid, standardize, control) {
    dX <- dim(X)

    scale.x <- 1
    mux <- 0
    muy <- 0

    if (standardize == TRUE) {
        scale.x <- apply(X, 2, mad)
        mux <- apply(X, 2, median)
        muy <- median(y)

        X <- scale(X, center = mux, scale = scale.x)
        y <- y - muy
    }

    # final.coefficients <- matrix(NA_real_, nrow = dX[2L] + 1L, ncol = length(lambda.grid))
    final.estimates <- vector("list", length(lambda.grid))

    ## For the first value of lambda, we will do a cold start
    init.current <- initest.cold(X, y, alpha, lambda.grid[1L], control)$initCoef

    for (i in seq_along(lambda.grid)) {
        full.all <- apply(init.current, 2, function(coef, lambda) {
            pen.s.reg(X, y, alpha, lambda,
                      maxit = control$pense.maxit,
                      init.coef = coef,
                      control)
        }, lambda = lambda.grid[i])

        best.est <- which.min(sapply(full.all, function(full) {
            full$objF
        }))
        init.current <- matrix(c(full.all[[best.est]]$intercept, full.all[[best.est]]$beta), ncol = 1L)

        final.estimates[[i]] <- full.all[[best.est]]

        if (standardize == TRUE) {
            final.estimates[[i]]$beta <- final.estimates[[i]]$beta / scale.x
            final.estimates[[i]]$intercept <- final.estimates[[i]]$intercept + muy -
                drop(final.estimates[[i]]$beta %*% mux)
        }
    }

    return(final.estimates)
}


initest.cold <- function(X, y, alpha, lambda, control) {
    lambda1 <- alpha * lambda
    lambda2 <- 0.5 * lambda * (1 - alpha)

    initraw <- enpy(X, y, lambda1, lambda2,
                    deltaesc = control$mscale.delta,
                    cc.scale = control$mscale.cc,
                    psc.method = control$init.psc.method,
                    prosac = control$init.psc.proportion,
                    clean.method = control$init.resid.clean.method,
                    C.res = control$init.resid.threshold,
                    prop = control$init.resid.proportion,
                    py.nit = control$init.maxit,
                    en.tol = control$init.tol,
                    control = enpy.control(
                        en.maxit = control$init.en.maxit,
                        en.tol = control$init.en.tol,
                        en.centering = TRUE,
                        mscale.maxit = control$mscale.maxit,
                        mscale.tol = control$mscale.tol,
                        mscale.rho.fun = control$mscale.rho.fun
                    ))

    ## Compute a "short" PENSE for each candidate solution
    initconc <- apply(initraw$coeff, 2, function(coef) {
        conc <- pen.s.reg(X, y, alpha, lambda,
                          maxit = control$init.csteps,
                          init.coef = coef,
                          control,
                          warn = FALSE)

        return(c(
            conc$objF,
            conc$intercept,
            conc$beta
        ))
    })

    objFRanking <- orderOmitTies(initconc[1L, ], tol = control$init.tol)

    nkeep <- min(length(objFRanking$index), control$init.nkeep)
    indkeep <- objFRanking$index[seq_len(nkeep)]

    return(list(
        initCoef = initconc[-1L, indkeep],
        objF = initconc[1L, indkeep]
    ))
}


