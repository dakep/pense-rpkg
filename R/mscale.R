#' Robust M-estimate of Scale
#'
#' Compute the M-estimate of scale with MAD as initial estimate.
#'
#' This solves the M-estimation equation given by
#' \deqn{\sum_{i=1}^n \rho( x_i / s_n; cc ) = n \delta}
#'
#' @note \code{NA} values in \code{x} are removed before calculating the
#'      M-scale.
#'
#' @param x numeric vector of observations.
#' @param delta target value of the M-estimation equation.
#' @param rho rho function to use in the M-estimation equation.
#' @param cc non-negative constant for the chosen rho function. If missing,
#'      it will be chosen such that the expected value of the rho function
#'      under the normal model is equal to \code{delta}.
#' @param eps threshold for convergence.
#' @param maxit maximum number of iterations.
#'
#' @return the M-scale as a numeric vector of length one.
#'
#' @example examples/mscale.R
#' @useDynLib pense, .registration = TRUE
#' @export
mscale <- function(x, delta = 0.5, rho = c("bisquare", "huber"), cc,
                   eps = 1e-8, maxit = 200) {

    if (!is.numeric(x) || !is.null(dim(x)) || length(x) == 0) {
        stop("`x` must be a numeric vector.")
    }
    if (any(is.na(x))) {
        x <- x[!is.na(x)]
    }

    rho <- .rho2IntRho(match.arg(rho))
    delta <- .check_arg(delta, "numeric", range = c(0, 1))

    cc <- if (missing(cc)) {
        consistency.rho(delta, rho)
    } else {
        .check_arg(cc, "numeric", range = 0)
    }

    eps <- .check_arg(eps, "numeric", range = 0)
    maxit <- .check_arg(maxit, "integer", range = 0)

    x <- as.numeric(x)

    scale <- .Call(
        C_mscale,
        x,
        length(x),
        delta,
        cc,
        maxit,
        eps,
        rho
    )

    if (!is.finite(scale)) {
        scale <- NA_real_
    }

    return(scale)
}
