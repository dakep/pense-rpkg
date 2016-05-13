##
#' @useDynLib penseinit C_augtrans C_pen_mstep
#' @importFrom robustbase Mchi
pensemstep <- function(X, y, cc, init.scale, init.coef, alpha, lambda, control) {
    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X, dX[1L], dX[2L])
    dX[2L] <- dX[2L] + 1L

    cctrl <- penseinit:::initest.control(
        numIt = control$pense.maxit,
        eps = control$pense.tol^2,
        mscale.delta = control$mscale.delta,
        mscale.cc = cc,
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

    res <- .Call(C_pen_mstep, Xtr, y, dX[1L], dX[2L], init.coef, init.scale, alpha, lambda, cctrl)

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

    mrhoinf <- 1 / MrhoInf(cc, control$mscale.rho.fun)

    repeat {
        it <- it + 1L

        resid.scaled <- resid / init.scale
        # Mwgt is safer then Mchi in the case 0/0!
        wbeta <- Mwgt(resid.scaled, cc = cc, psi = control$mscale.rho.fun) * mrhoinf

        Wbeta.tilde <- sqrt(wbeta)

        Xweight <- X * Wbeta.tilde
        yweight <- (y - current.coefs[1L]) * Wbeta.tilde

        beta.obj <- .elnet.fit(Xweight, yweight, alpha = alpha, lambda = lambda,
                               centering = FALSE, maxit = control$pense.en.maxit,
                               eps = control$pense.en.tol, warmCoef = current.coefs,
                               en.algorithm = control$en.algorithm)

        prev.coefs <- current.coefs
        current.coefs <- beta.obj$coefficients

        resid <- as.vector(y - X %*% current.coefs[-1L])
        # New intercept
        wintercept <- wbeta / sum(wbeta)
        current.coefs[1L] <- sum(wintercept * resid)

        # Update residuals and scale
        resid <- resid - current.coefs[1L]
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


################################################################################
##
## All of the following code is taken from R package mmlasso-1.3.4
## Copyright Ezequiel Smucler <ezequiels.90@gmail.com>
##
##################################################################################
#' @importFrom mmlasso mmlasso
pensemstepL1 <- function(Xs, y, cc, init.scale, init.coef, alpha, lambda, control) {
    #Performs iteratively re-weighted Lasso
    #INPUT
    #Xs,y: data
    #init.scale, init.coef: initial estimate of regression and scale
    #lambda: penalization parameter
    #c1: tuning constant for the rho function
    #niter.mm: maximum number of iterations
    #OUTPUT:
    #coef: final estimate

    MMcpp_ini <- MMLassoCpp_ini(Xs, y, init.coef) # Compute residuals and augment Xs with column of 1's
    x <- MMcpp_ini$x # Xs with a column of 1's in front
    resid.n <- MMcpp_ini$resid.n # Residuals for init.coef
    tol <- 1
    m <- 0L
    beta.n <- init.coef
    p <- length(init.coef)
    u.Gram <- !(p >= 500)

    while (tol >= control$pense.tol){
        beta.o <- beta.n
        if(all(beta.o == 0)){
            return(beta.o)
        }
        MMcpp1 <- MMLassoCpp1(x, y, beta.o, init.scale, cc)
        xort <- MMcpp1$xort # X* orthogonal in the paper
        xjota <- MMcpp1$xjota # vector with sqrt(weight_i)
        yast <- MMcpp1$yast # y* = W %*% y, W = diag(sqrt(weight_i))
        alpha <- MMcpp1$alpha # .. eta_j in the paper
        if(all(xjota == 0)){
            return(beta.o)
        }

        res.elnet <- .elnet.fit(xort, yast, alpha = 1, lambda = lambda,
                                maxit = control$pense.en.maxit, eps = control$pense.en.tol,
                                centering = FALSE, warmCoef = beta.o,
                                en.algorithm = control$en.algorithm)
        beta.Lasso <- res.elnet$coefficients[-1L]

        # try(res.Lasso <- robustHD:::fastLasso(xort, yast, 2 * lambda, intercept=FALSE,
        #                                       normalize=FALSE, use.Gram = u.Gram), silent = TRUE)
        # if (class(res.Lasso)=="try-error"){
        #     warning("fastLasso failed")
        #     return(beta.o)
        # }
        # beta.Lasso <- res.Lasso$coefficients


        MMcpp2 <- MMLassoCpp2(xjota, yast, beta.Lasso, beta.o, alpha) # compute intercept and rel. change
        beta.n <- MMcpp2$beta.n # beta.Lasso + intercept
        tol <- MMcpp2$tol # rel. change from previous solution
        m <- m + 1L
        if (m >= control$pense.maxit) {
            break
        }
    }


    ret <- list(
        intercept = beta.n[1L],
        beta = beta.n[-1L],
        resid = NULL,
        rel.change = tol,
        iterations = m,
        objF = NA_real_
    )


    ret$resid <- drop(y - ret$intercept - Xs %*% ret$beta)

    ret$objF <- sum(Mchi(drop(ret$resid / init.scale), cc = cc, psi = control$mscale.rho.fun)) +
        lambda * sum(abs(ret$beta))

    if (ret$rel.change > control$pense.tol) {
        warning(sprintf("M-step did not converge for lambda = %.3f", lambda))
    }

    return(ret)
}


MMLassoCpp_ini <- function(xx, y, beta_ini) {
    .Call('mmlasso_MMLassoCpp_ini', PACKAGE = 'mmlasso', xx, y, beta_ini)
}

MMLassoCpp1 <- function(x, y, beta_ini, scale_ini, c1) {
    .Call('mmlasso_MMLassoCpp1', PACKAGE = 'mmlasso', x, y, beta_ini, scale_ini, c1)
}

MMLassoCpp2 <- function(xjota, yast, beta_lars, beta_o, alpha) {
    .Call('mmlasso_MMLassoCpp2', PACKAGE = 'mmlasso', xjota, yast, beta_lars, beta_o, alpha)
}
