#' PY (Peña-Yohai) initial estimates for EN
#'
#' Computes the PY initial estimates for EN with different strategies for computing
#' the principal sensitivity components
#'
#' Three different methods to calculate the sensitivity components are implemented:
#'
#' @param X The data matrix X
#' @param y The response vector
#' @param lambda1,lambda2 The EN penalty parameters (adjusted for the number of observations
#'          in \code{X})
#' @param deltasc,cc.scale Parameters for the M-equation of the scale. Tukey's bisquare rho function
#'      is used internally.
#' @param psc.method The method to use for computing the principal sensitivity components.
#'      See details.
#' @param prosac The proportion of observations to remove based on PSCs
#' @param clean.method How to clean the data based on large residuals.
#'          If \code{"threshold"}, all observations with scaled residuals larger than
#'          \code{C.res} will be removed.
#'          If \code{"proportion"}, observations with the largest \code{prop} residuals
#'          will be removed.
#' @param py.nit The maximum number of iterations to perform.
#' @param en.tol The relative tolerance for convergence.
#'
#' @return \item{initCoef}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
#'
#' @export
enpy <- function(X, y, lambda1, lambda2, deltasc, cc.scale,
                 psc.method=c("rr","Qp","Mn"), prosac,
                 clean.method = c("threshold", "proportion"), C.res, prop,
                 py.nit, en.tol) {
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

    psc.method <- match.arg(psc.method)

    ## Add leading column of 1's
    X <- cbind(1, X)

    switch(psc.method,
           rr = enpy.rr(X, y, lambda1, lambda2, deltasc, cc.scale, prosac, clean.method, C.res,
                        prop, py.nit, en.tol),
           Qp = stop("Method `Qp` is not yet implemented"),
           Mn = stop("Method `Mn` is not yet implemented"))
}
