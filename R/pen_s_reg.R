## Internal function to calculate PENSE for a given initial estimate (init.coef)
## in C++
##
#' @useDynLib pense C_augtrans C_pen_s_reg
pen_s_reg <- function(X, y, alpha, lambda, init_coef, warn = TRUE,
                      options, en_options) {
    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X)
    dX[2L] <- dX[2L] + 1L

    lambda1 <- alpha * lambda
    lambda2 <- lambda * (1 - alpha) / 2

    res <- .Call(
        C_pen_s_reg,
        Xtr,
        y,
        init_coef,
        alpha,
        lambda,
        options,
        en_options
    )

    ret <- list(
        intercept = res[[1L]][1L],
        beta = res[[1L]][-1L],
        resid = res[[2L]],
        scale = res[[3L]],
        rel_change = sqrt(res[[4L]]),
        iterations = res[[5L]],
        objF = NA_real_
    )

    ret$objF <- ret$scale^2 + lambda * (
        0.5 * (1 - alpha) * sum(ret$beta^2) + alpha * sum(abs(ret$beta))
    )

    ##
    ## Check if the S-step converged.
    ## Be extra careful with the comparison as rel_change may be NaN or NA
    ##
    if (!isTRUE(ret$rel_change < options$eps) && isTRUE(warn)) {
        warning(sprintf("PENSE did not converge for lambda = %.6f", lambda))
    }

    return(ret)
}



## Internal function to calculate PENSE for a given initial estimate
## (init.coef) in R
##
#' @importFrom robustbase Mwgt MrhoInf
pen_s_reg_rimpl <- function(X, y, alpha, lambda, init_coef, warn = TRUE,
                            options, en_options) {
    dX <- dim(X)
    p <- dX[2L]
    n <- dX[1L]

    current.coefs <- init_coef
    prev.coefs <- NULL

    resid <- as.vector(y - current.coefs[1L] - X %*% current.coefs[-1L])
    scale <- mscale(
        resid,
        cc = options$cc,
        b = options$bdp,
        eps = control$mscaleEps,
        maxit = control$mscaleMaxit
    )

    tol <- control$eps^2
    it <- 0L
    rel_change <- Inf


    repeat {
        it <- it + 1L

        resid.scaled <- resid / scale
        # Mwgt is safer then Mchi in the case 0/0!
        wbeta <- Mwgt(resid.scaled, cc = control$cc, psi = "bisquare")
        # Not necessary to compute the "true" weights, as the normalization
        # will remove this again
        # weights <- wbeta * (scale / sum(resid^2 * wbeta))

        # normalize weights to sum to n
        weights <- n * wbeta / sum(wbeta)

        ## Perform weighted elastic net
        weight.en <- .elnet.wfit(
            X,
            y,
            weights = weights,
            alpha = alpha,
            lambda = 0.5 * lambda,
            en_options
        )

        prev.coefs <- current.coefs
        current.coefs <- weight.en$coefficients

        resid <- weight.en$residuals

        scale <- mscale(
            resid,
            cc = options$cc,
            b = options$bdp,
            eps = control$mscaleEps,
            maxit = control$mscaleMaxit
        )

        rel_change <- sum((prev.coefs - current.coefs)^2) / sum(prev.coefs^2)

        if (is.nan(rel_change)) { # if 0/0
            rel_change <- 0
        }

        ## Check convergence
        if (it >= maxit || all(rel_change < tol)) {
            break
        }
    }

    ##
    ## Check if the S-step converged.
    ## Be extra careful with the comparison as rel.change may be NaN or NA
    ##
    if (!isTRUE(rel_change < tol) && isTRUE(warn)) {
        warning(sprintf("PENSE did not converge for lambda = %.6f", lambda))
    }

    beta <- current.coefs[-1L]
    objf <- scale^2 + lambda * (
        0.5 * (1 - alpha) * sum(beta^2) + alpha * sum(abs(beta))
    )

    return(list(
        intercept = current.coefs[1L],
        beta = beta,
        resid = resid,
        scale = scale,
        iterations = it,
        rel_change = sqrt(rel_change),
        objF = objf
    ))
}
