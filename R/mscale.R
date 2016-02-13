#' Robust M-estimate of Scale
#'
#' Compute the M-estimate of scale with the MAD as initial estimate.
#'
#' This solves the M-estimation equation given by
#' \deqn{\sum_{i=1}^n \rho( x_i / s_n; cc ) = n b}
#'
#' \code{NA} values in \code{x} are removed before calculating the
#'
#' @param x A numeric vector.
#' @param b The value of the M-estimation equation.
#' @param rho The rho function to use in the M-estimation equation.
#' @param cc Non-negative constant for the chosen rho function.
#' @param eps Threshold for convergence
#' @param max.it The maximum number of iterations
#'
#' @return Numeric vector of length one
#'
#' @useDynLib penseinit
#' @export
mscale <- function(x, b = 0.5, rho = c("bisquare", "huber", "gauss"), cc = 1.54764,
                   eps = 1e-8, max.it = 200) {

    if (!is.numeric(x) || !is.null(dim(x)) || length(x) == 0) {
        stop("`x` must be a numeric vector.")
    }
    if (any(is.na(x))) {
        x <- x[!is.na(x)]
    }

    ctrl <- initest.control(mscaleB = b,
                            mscaleCC = cc,
                            mscaleMaxIt = max.it,
                            mscaleEPS = eps,
                            mscaleRhoFun = rho,
                            ## Not needed arguments below
                            lambda1 = 0,
                            lambda2 = 0,
                            numIt = 1,
                            residThreshold = 1,
                            residProportion = 0.5,
                            pscProportion = 0.5)

    x <- as.numeric(x)

    scale <- .Call("C_mscale", x, length(x), ctrl$mscaleB, ctrl$mscaleCC, ctrl$mscaleMaxIt,
                   ctrl$mscaleEPS, ctrl$mscaleRhoFun, PACKAGE = "penseinit")

    if (!is.finite(scale)) {
        scale <- NA_real_
    }

    return(scale)
}
