#' PY (Pena-Yohai) initial estimates for EN
#'
#' Computes the PY initial estimates for EN with different strategies for computing
#' the principal sensitivity components
#'
#' Three different methods to calculate the sensitivity components are implemented:
#' \describe{
#'      \item{\code{"rr"}}{Approximate the PSCs by using the residuals from the elastic net fit
#'                         and the hat matrix from the ridge regression. This method only works
#'                         if \code{alpha} < 1 or \code{ncol(X)} < \code{nrow(X)}.}
#'      \item{\code{"Mn"}}{Calculate the PSCs from the difference between the
#'                         residuals and leave-one-out residuals from elastic net.}
#' }
#'
#' @param X The data matrix X.
#' @param y The response vector.
#' @param alpha,lambda The EN penalty parameters (NOT adjusted for the number of observations
#'          in \code{X}).
#' @param deltaesc,cc.scale Parameters for the M-equation of the scale. If \code{cc.scale}
#'          is missing or invalid, it will be chosen such that the expected value of the
#'          rho function under the normal model is equal to \code{deltaesc}.
#'          The default rho function (Tukey's bisquare) can be changed by
#'          parameter \code{control}.
#' @param psc.method The method to use for computing the principal sensitivity components.
#'      See details.
#' @param prosac The proportion of observations to remove based on PSCs.
#' @param clean.method How to clean the data based on large residuals.
#'          If \code{"threshold"}, all observations with scaled residuals larger than
#'          \code{C.res} will be removed.
#'          If \code{"proportion"}, observations with the largest \code{prop} residuals
#'          will be removed.
#' @param C.res,prop See \code{clean.method} for details.
#' @param py.nit The maximum number of iterations to perform.
#' @param en.tol The relative tolerance for convergence.
#' @param control Optional further control parameters from \code{\link{enpy.control}}.
#'
#' @return \item{initCoef}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
#'
#' @export
enpy <- function(X, y, alpha, lambda, deltaesc, cc.scale,
                 psc.method=c("rr", "Mn", "Qp"), prosac,
                 clean.method = c("threshold", "proportion"), C.res = NULL, prop = NULL,
                 py.nit, en.tol, control = enpy.control()) {
    y <- drop(y)

    dX <- dim(X)
    dY <- dim(y)
    yl <- length(y)

    if (is.null(yl) || (!is.null(dY) && length(dY) != 1L) || !is.numeric(y)) {
        stop("`yl` must be a numeric vector")
    }

    if (is.null(dX) || length(dX) != 2L || !is.numeric(X) || dX[1L] != yl) {
        stop("`X` must be a numeric matrix with the same number of observations as `y`")
    }

    if (anyNA(X) || anyNA(y)) {
        stop("Missing values are not supported")
    }

    if (missing(cc.scale) || length(cc.scale) != 1L || !is.numeric(cc.scale) || is.na(cc.scale) ||
        cc.scale <= 0) {
        cc.scale <- 0
    }

    psc.method <- match.arg(psc.method)

    if (alpha == 0) {
        psc.method <- 'rr'
    }

    result <- switch(psc.method,
                     rr = enpy.rr(X, y, alpha, lambda, deltaesc, cc.scale, prosac, clean.method,
                                  C.res, prop, py.nit, en.tol, control),
                     Mn = enpy.exact(X, y, alpha, lambda, deltaesc, cc.scale,
                                     psc.method, prosac, clean.method, C.res, prop, py.nit,
                                     en.tol, control),
                     stop(sprintf("`psc.method` '%s' not supported", psc.method)))

    resorder <- sort.list(result$objF, na.last = NA, method = "quick")

    dups <- which(diff(result$objF[resorder]) < en.tol)
    if (length(dups) > 0) {
        remove <- resorder[dups]
        result$objF <- result$objF[-remove]
        result$coeff <- result$coeff[ , -remove, drop = FALSE]
    }

    return(result)
}
