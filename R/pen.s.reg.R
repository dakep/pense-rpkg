## Internal function to calculate PENSE for a given initial estimate (init.coef)
## in C++
##
#' @useDynLib penseinit C_augtrans C_pen_s_reg
pen.s.reg <- function(X, y, alpha, lambda, init.coef, maxit, control, warn = TRUE) {
    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X, dX[1L], dX[2L])
    dX[2L] <- dX[2L] + 1L

    lambda1 <- alpha * lambda * dX[1L]
    lambda2 <- lambda * (1 - alpha) * dX[1L] / 2

    cctrl <- penseinit:::initest.control(
        numIt = maxit,
        eps = control$pense.tol^2,
        mscale.delta = control$mscale.delta,
        mscale.cc = control$mscale.cc,
        enpy.control = enpy.control(
            en.maxit = control$pense.en.maxit,
            en.tol = control$pense.en.tol,
            en.centering = FALSE,
            mscale.maxit = control$mscale.maxit,
            mscale.tol = control$mscale.tol,
            mscale.rho.fun = "bisquare",
            en.algorithm = control$en.algorithm
        ),

        # Not needed, set for completeness
        lambda1 = 0,
        lambda2 = 0,
        resid.clean.method = "proportion",
        resid.threshold = 0.5,
        resid.proportion = 0.5,
        psc.proportion = 0.5
    )

    res <- .Call(C_pen_s_reg, Xtr, y, dX[1L], dX[2L], init.coef, alpha, lambda, cctrl)

    ret <- list(
        intercept = res[[1L]][1L],
        beta = res[[1L]][-1L],
        resid = res[[2L]],
        scale = res[[3L]],
        rel.change = sqrt(res[[4L]]),
        iterations = res[[5L]],
        objF = NA_real_
    )

    ret$objF <- dX[1L] * ret$scale^2 + lambda1 * sum(abs(ret$beta)) + lambda2 * sum(ret$beta^2)

    if (ret$rel.change > control$pense.tol && identical(warn, TRUE)) {
        warning(sprintf("PENSE did not converge for lambda = %.3f", lambda))
    }

    return(ret)
}



## Internal function to calculate PENSE for a given initial estimate (init.coef) in R
##
#' @importFrom robustbase Mwgt MrhoInf
pen.s.reg.rimpl <- function(X, y, alpha, lambda, init.coef, maxit, control, warn = TRUE) {
    dX <- dim(X)
    p <- dX[2L]
    n <- dX[1L]

    current.coefs <- init.coef
    prev.coefs <- NULL

    resid <- as.vector(y - current.coefs[1L] - X %*% current.coefs[-1L])
    scale <- mscale(resid, cc = control$mscale.cc,
                    b = control$mscale.delta, rho = control$mscale.rho.fun,
                    eps = control$mscale.tol, max.it = control$mscale.maxit)

    lambda1 <- alpha * lambda * n
    lambda2 <- lambda * (1 - alpha) * n / 2

    tol <- control$pense.tol^2
    it <- 0L
    rel.change <- Inf

    mrhoinf <- 1 / MrhoInf(control$mscale.cc, control$mscale.rho.fun)
    en.eps <- control$pense.en.tol

    repeat {
        it <- it + 1L

        resid.scaled <- resid / scale
        # Mwgt is safer then Mchi in the case 0/0!
        wbeta <- Mwgt(resid.scaled, cc = control$mscale.cc, psi = control$mscale.rho.fun) * mrhoinf
        tau.beta <- scale^2 / sum(resid^2 * wbeta)
        Wbeta.tilde <- sqrt(tau.beta * wbeta)

        Xweight <- X * Wbeta.tilde
        yweight <- (y - current.coefs[1L]) * Wbeta.tilde

        ## Adjust convergence threshold for elastic net if we use coordinate descent
        if (control$en.algorithm == "coordinate-descent") {
            en.eps <- control$pense.en.tol * dX[1L] * mscale(yweight)^2
        }

        ## We need to
        beta.obj <- .elnet.fit(Xweight, yweight, alpha = alpha, lambda = lambda / (2 * n),
                               centering = FALSE, maxit = control$pense.en.maxit,
                               eps = en.eps, warmCoef = current.coefs,
                               en.algorithm = control$en.algorithm)

        prev.coefs <- current.coefs
        current.coefs <- beta.obj$coefficients

        resid <- as.vector(y - X %*% current.coefs[-1L])
        # New intercept
        wintercept <- wbeta / sum(wbeta)
        current.coefs[1L] <- sum(wintercept * resid)

        # Update residuals and scale
        resid <- resid - current.coefs[1L]
        scale <- mscale(resid, cc = control$mscale.cc,
                        b = control$mscale.delta, rho = control$mscale.rho.fun,
                        eps = control$mscale.tol, max.it = control$mscale.maxit)

        rel.change <- sum((prev.coefs - current.coefs)^2) / sum(prev.coefs^2)

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

    beta <- current.coefs[-1L]
    objf <- n * scale^2 + lambda2 * sum(beta^2) + lambda1 * sum(abs(beta))

    return(list(
        intercept = current.coefs[1L],
        beta = beta,
        resid = resid,
        scale = scale,
        iterations = it,
        rel.change = sqrt(rel.change),
        objF = objf
    ))
}
