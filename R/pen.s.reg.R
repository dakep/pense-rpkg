## Internal function to calculate PENSE for a given initial estimate (init.coef)
## @param mscale.cc is used instead of control$mscale.cc
##
#' @importFrom robustbase Mchi
pen.s.reg <- function(X, y, maxit, lambda1, lambda2, init.coef, mscale.cc, control) {
    dX <- dim(X)
    p <- dX[2L]
    n <- dX[1L]

    #Intercept, slopes, and scale
    beta <- init.coef[-1L]
    intercept <- init.coef[1L]

    prev.beta <- NULL

    resid <- as.vector(y - intercept - X %*% beta)
    scale <- mscale(resid, cc = mscale.cc,
                    b = control$mscale.delta, rho = control$mscale.rho.fun,
                    eps = control$mscale.tol, max.it = control$mscale.maxit)

    enlambda <- (2 * lambda2 + lambda1) / n
    enalpha <- lambda1 / (2 * lambda2 + lambda1)
    if (!is.finite(enalpha)) {
        enalpha <- 0
    }

    tol <- control$pense.tol^2
    it <- 0L
    rel.change <- Inf

    repeat {
        it <- it + 1L

        resid.scaled <- resid / scale
        wbeta <- Mchi(resid.scaled, cc = mscale.cc, deriv = 1,
                      psi = control$mscale.rho.fun) / resid.scaled
        tau.beta <- n * 2 * scale^2 / sum(resid^2 * wbeta)
        Wbeta.tilde <- sqrt(tau.beta * wbeta)

        Xweight <- X * Wbeta.tilde
        yweight <- (y - intercept) * Wbeta.tilde

        beta.obj <- .elnet.fit(Xweight, yweight, alpha = enalpha, lambda = enlambda,
                               centering = FALSE, maxit = control$pense.en.maxit,
                               eps = control$pense.en.tol)

        prev.beta <- beta
        beta <- beta.obj$coefficients[-1L]

        resid <- as.vector(y - X %*% beta)
        # New intercept
        wintercept <- wbeta / sum(wbeta)
        intercept <- sum(wintercept * resid)

        # Update residuals and scale
        resid <- resid - intercept
        scale <- mscale(resid, cc = mscale.cc,
                        b = control$mscale.delta, rho = control$mscale.rho.fun,
                        eps = control$mscale.tol, max.it = control$mscale.maxit)

        rel.change <- sum((prev.beta - beta)^2) / sum(prev.beta^2)

        ## Check convergence
        if (it >= maxit || rel.change < tol) {
            break
        }
    }

    if (rel.change > tol) {
        warning("PENSE did not converge")
    }

    objf <- n * scale^2 + lambda2 * sum(beta^2) + lambda1 * sum(abs(beta))

    return(list(
        intercept = intercept,
        beta = beta,
        resid = resid,
        scale = scale,
        iterations = it,
        rel.change = rel.change,
        objF = objf
    ))
}
