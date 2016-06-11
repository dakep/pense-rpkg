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
