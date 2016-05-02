################################################################################
##
## All of the following code is taken from R package mmlasso-1.3.4
## Copyright Ezequiel Smucler <ezequiels.90@gmail.com>
##
##################################################################################


# @-import robustHD
#' @importFrom mmlasso mmlasso
pensemstep <- function(Xs, y, cc, init.scale, init.coef, lambda, control) {
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
    return(beta.n)
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
