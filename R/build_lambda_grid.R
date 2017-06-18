## Generate a grid of lambda values
## The grid is *INDEPENDENT* from the sample size
##
#' @importFrom stats mad median
#' @importFrom robustbase scaleTau2 covGK
build_lambda_grid <- function(X, y, alpha, nlambda, standardize = TRUE,
                              lambda_min_ratio = NULL, scaleFun = mad,
                              meanFun = median) {
    dX <- dim(X)

    if (isTRUE(standardize)) {
        scale.x <- apply(X, 2, scaleFun)
        mu.x <- apply(X, 2, meanFun)

        Xs <- scale(X, scale = scale.x, center = mu.x)
        y <- y - meanFun(y)
    } else {
        scale.x <- rep.int(1, dX[2L])
        Xs <- X
    }

    # Compute a grid of lambdas for the EN estimator
    ycs <- y / scaleFun(y)

    ## Pairwise robust covariances
    covxy <- apply(Xs, 2, covGK, ycs, sigmamu = scaleTau2)

    if (is.null(lambda_min_ratio)) {
        lambda_min_ratio <- min(1e-4, 1e-4 * 10^floor(log10(dX[2L] / dX[1L])))
    }

    lmax <- max(abs(covxy))
    lmin <- lambda_min_ratio * lmax

    lambda <- exp(seq(log(lmin), log(2 * lmax), length.out = nlambda)) *
        2 * scaleFun(y) / alpha
    #adjustment for an S-est (may not be enough)
    # lambda <- lambda * 1.1 * 2

    return(lambda)
}
