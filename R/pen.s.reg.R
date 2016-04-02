## Internal function to calculate PENSE for a given initial estimate (init.coef)
##
#' @importFrom robustbase Mchi
pen.s.reg <- function(X, y, alpha, lambda, init.coef, maxit, control, warn = TRUE) {
    dX <- dim(X)
    p <- dX[2L]
    n <- dX[1L]

    #Intercept, slopes, and scale
    beta <- init.coef[-1L]
    intercept <- init.coef[1L]

    prev.beta <- NULL

    resid <- as.vector(y - intercept - X %*% beta)
    scale <- mscale(resid, cc = control$mscale.cc,
                    b = control$mscale.delta, rho = control$mscale.rho.fun,
                    eps = control$mscale.tol, max.it = control$mscale.maxit)

    lambda1 <- alpha * lambda * n
    lambda2 <- lambda * (1 - alpha) * n / 2

    tol <- control$pense.tol^2
    it <- 0L
    rel.change <- Inf

    repeat {
        it <- it + 1L

        resid.scaled <- resid / scale
        wbeta <- Mchi(resid.scaled, cc = control$mscale.cc, deriv = 1,
                      psi = control$mscale.rho.fun) / resid.scaled
        tau.beta <- n * 2 * scale^2 / sum(resid^2 * wbeta)
        Wbeta.tilde <- sqrt(tau.beta * wbeta)

        Xweight <- X * Wbeta.tilde
        yweight <- (y - intercept) * Wbeta.tilde

        ## Adjust convergence threshold for elastic net
        en.eps <- control$pense.en.tol * dX[1L] * mscale(yweight)^2

        beta.obj <- .elnet.fit(Xweight, yweight, alpha = alpha, lambda = lambda,
                               centering = FALSE, maxit = control$pense.en.maxit,
                               eps = control$pense.en.tol)

        prev.coefs <- c(intercept, beta)
        beta <- beta.obj$coefficients[-1L]

        resid <- as.vector(y - X %*% beta)
        # New intercept
        wintercept <- wbeta / sum(wbeta)
        intercept <- sum(wintercept * resid)

        # Update residuals and scale
        resid <- resid - intercept
        scale <- mscale(resid, cc = control$mscale.cc,
                        b = control$mscale.delta, rho = control$mscale.rho.fun,
                        eps = control$mscale.tol, max.it = control$mscale.maxit)

        rel.change <- sum((prev.coefs - c(intercept, beta))^2) / sum(prev.coefs^2)

        if (is.nan(rel.change)) { # if 0/0
            rel.change <- 0
        }

        ## Check convergence
        if (it >= maxit || all(rel.change < tol)) {
            break
        }
    }

    if (rel.change > tol && identical(warn, TRUE)) {
        warning(sprintf("PENSE did not converge for lambda = %.3f", lambda))
    }

    objf <- n * scale^2 + lambda2 * sum(beta^2) + lambda1 * sum(abs(beta))

    return(list(
        intercept = intercept,
        beta = beta,
        resid = resid,
        scale = scale,
        iterations = it,
        rel.change = sqrt(rel.change),
        objF = objf
    ))
}
