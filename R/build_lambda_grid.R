## Generate a grid of lambda values
## The grid is *INDEPENDENT* from the sample size
##
#' @importFrom stats mad median
#' @importFrom robustbase scaleTau2 covGK
build_lambda_grid <- function(x, y, alpha, nlambda, lambda_min_ratio) {
    dx <- dim(x)

    max_scale_x <- max(apply(x, 2, mad))

    # Compute a grid of lambdas for the EN estimator
    scale_y <- mad(y)
    ys <- y / scale_y

    ## Pairwise robust covariances
    covxy <- apply(x, 2, covGK, ys, sigmamu = scaleTau2)

    lmax <- max(abs(covxy)) * 2 * scale_y / (max_scale_x * max(0.01, alpha))
    lmin <- lambda_min_ratio * lmax

    lambda <- exp(seq(log(lmin), log(lmax), length.out = nlambda))

    return(lambda)
}
