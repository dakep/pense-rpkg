#' @importFrom robustbase Mchi
pen.s.reg <- function(x, y, maxit, lambda1, lambda2, init.coef, control) {
    dX <- dim(X)
    p <- dX[2L]
    n <- dX[1L]

    Rhop <- function(r) {
        Mchi(r, cc = control$mscale.cc, deriv = 1, psi = control$mscale.rho.fun)
    }

    #Intercept, slopes, and scale
    beta <- init.coef[-1L]
    intercept <- init.coef[1L]

    prev.beta <- NULL

    resid <- as.vector(y - intercept - X %*% beta)
    scale <- mscale(resid, b = control$mscale.delta, rho = control$mscale.rho.fun,
                    cc = control$mscale.cc, eps = control$mscale.tol,
                    max.it = control$mscale.maxit)

    enlambda <- (2 * lambda2 + lambda1) / n
    enalpha <- lambda1 / (2 * lambda2 + lambda1)

    tol <- control$pense.en.tol^2
    it <- 1L

    repeat {
        resid.scaled <- resid / scale
        wbeta <- Mchi(resid.scaled, cc = control$mscale.cc, deriv = 1,
                      psi = control$mscale.rho.fun) / resid.scaled
        tau.beta <- n * 2 * scale^2 / sum(resid^2 * wbeta)
        Wbeta.tilde <- sqrt(tau.beta * wbeta)

        Xweight <- X * Wbeta.tilde
        yweight <- (y - intercept) * Wbeta.tilde

        beta.obj <- elnet(Xweight, yweight, alpha = enalpha, lambda = enlambda,
                          centering = FALSE, maxit = control$pense.en.maxit)

        prev.beta <- beta
        beta <- beta.obj$coefficients[-1L]

        resid <- as.vector(y - X %*% beta)
        # New intercept
        wintercept <- wbeta / sum(wbeta)

        intercept <- sum(wintercept * resid)

        # Update residuals and scale
        resid <- resid - intercept
        scale <- mscale(resid, b = control$mscale.delta, rho = control$mscale.rho.fun,
                        cc = control$mscale.cc, eps = control$mscale.tol,
                        max.it = control$mscale.maxit)

        rel.change <- sum((prev.beta - beta)^2) / sum(prev.beta^2)

        ## Check convergence
        if (it > maxit || rel.change > tol) {
            break
        }
    }

    objf <- n * scale^2 + lambda2 * sum(beta^2) + lambda1 * sum(abs(beta))

    return(list(
        intercept = intercept,
        beta = beta,
        resid = resid,
        scale = scale,
        iterations = it,
        objF = objf
    ))
}
