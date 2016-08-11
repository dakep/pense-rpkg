## Generate a grid of lambda values
## The grid is *INDEPENDENT* from the sample size
##
#' @importFrom stats mad median
#' @importFrom robustbase scaleTau2 covGK
lambda.grid <- function(X, y, alpha, nlambda, control, standardize = TRUE, lambda.min.ratio = NULL,
                        scaleFun = mad, meanFun = median) {
    dX <- dim(X)

    if (identical(standardize, TRUE)) {
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

    if (is.null(lambda.min.ratio)) {
        lambda.min.ratio <- min(1e-2, 1e-3 * 10^floor(log10(dX[2L] / dX[1L])))
    }

    lmax <- max(abs(covxy))
    lmin <- lambda.min.ratio * lmax

    ## the lambda grid is not yet adjusted for the sample size! This will be done separately
    ## during estimation!
    lambda <- exp(seq(log(lmin), log(lmax), length.out = nlambda)) * 2 * scaleFun(y) / alpha
    #adjustment for an S-est (may not be enough)
    lambda <- lambda * 1.1 * facon(control$mscale.delta) / control$mscale.cc^2

    return(lambda)
}
