#' Robust M-estimate of Scale
#'
#' Compute the M-estimate of scale with MAD as initial estimate.
#'
#' This solves the M-estimation equation given by
#' \deqn{\sum_{i=1}^n \rho( x_i / s_n; cc ) = n b}
#'
#' \code{NA} values in \code{x} are removed before calculating the
#'
#' @param x A numeric vector.
#' @param b The value of the M-estimation equation.
#' @param rho The rho function to use in the M-estimation equation.
#' @param cc Non-negative constant for the chosen rho function. If missing, it will be
#'          chosen such that the expected value of the rho function under the normal model
#'          is equal to \code{b}.
#' @param eps Threshold for convergence
#' @param max.it The maximum number of iterations
#'
#' @return Numeric vector of length one
#'
#' @useDynLib penseinit C_mscale
#' @export
mscale <- function(x, b = 0.5, rho = c("bisquare", "huber"), cc,
                   eps = 1e-8, max.it = 200) {

    if (!is.numeric(x) || !is.null(dim(x)) || length(x) == 0) {
        stop("`x` must be a numeric vector.")
    }
    if (any(is.na(x))) {
        x <- x[!is.na(x)]
    }

    if (missing(cc)) {
        cc <- consitency.rho(b, rho)
    }

    ctrl <- initest.control(mscale.delta = b,
                            mscale.cc = cc,
                            enpy.control = enpy.control(mscale.maxit = max.it,
                                                        mscale.tol = eps,
                                                        mscale.rho.fun = rho),
                            ## Not needed arguments below
                            lambda1 = 0,
                            lambda2 = 0,
                            numIt = 1,
                            resid.threshold = 1,
                            resid.proportion = 0.5,
                            psc.proportion = 0.5)

    x <- as.numeric(x)

    scale <- .Call(C_mscale, x, length(x), ctrl$mscale.delta, ctrl$mscale.cc, ctrl$mscale.maxit,
                   ctrl$mscale.tol, ctrl$mscale.rho.fun)

    if (!is.finite(scale)) {
        scale <- NA_real_
    }

    return(scale)
}
