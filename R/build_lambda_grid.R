## Generate a grid of lambda values
## The grid is *INDEPENDENT* from the sample size
##
#' @importFrom stats mad median
#' @importFrom robustbase scaleTau2 covGK
build_lambda_grid <- function(X, y, alpha, nlambda, lambda_min_ratio = NULL) {
    dX <- dim(X)

    max_scale_x <- max(apply(X, 2, mad))

    # Compute a grid of lambdas for the EN estimator
    scale_y <- mad(y)
    ys <- y / scale_y

    ## Pairwise robust covariances
    covxy <- apply(X, 2, covGK, ys, sigmamu = scaleTau2)

    if (is.null(lambda_min_ratio)) {
        lambda_min_ratio <- min(1e-5, 1e-5 * 10^floor(log10(dX[2L] / dX[1L])))
    }

    lmax <- max(abs(covxy)) * 2 * scale_y / (max_scale_x * alpha)
    lmin <- lambda_min_ratio * lmax

    lambda <- exp(seq(log(lmin), log(lmax), length.out = nlambda))
    #adjustment for an S-est (may not be enough)
    # lambda <- lambda * 1.1 * 2

    return(lambda)
}
