#' PY (Peña-Yohai) initial estimates for EN
#'
#' Computes the PY initial estimates for EN with exact principal sensitivity components
#'
#' @param X The data matrix X, with leading column of 1's.
#' @param y The response vector.
#' @param lambda1,lambda2 The EN penalty parameters (adjusted for the number of observations
#'          in \code{X}).
#' @param deltaesc,cc.scale Parameters for the M-equation of the scale. The default
#'          rho function is Tukey's bisquare. This can be changed by the parameter \code{control}.
#' @param psc.method The method to use for computing the principal sensitivity components.
#'          See details.
#' @param prosac The proportion of observations to remove based on PSCs.
#' @param clean.method How to clean the data based on large residuals.
#'          If \code{"threshold"}, all observations with scaled residuals larger than
#'          \code{C.res} will be removed.
#'          If \code{"proportion"}, observations with the largest \code{prop} residuals
#'          will be removed.
#' @param C.res,prop How many or which observations should be removed based on their residuals.
#' @param py.nit The maximum number of iterations to perform.
#' @param en.tol The relative tolerance for convergence.
#' @param control Further control parameters (see \code{\link{enpy.control}})
#'
#' @return \item{coeff}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
#'
#' @useDynLib penseinit C_enpy_exact
#' @importFrom Rcpp evalCpp
enpy.exact <- function(X, y, lambda1, lambda2, deltaesc, cc.scale,
                       psc.method = c("Mn", "Qp"), prosac,
                       clean.method = c("threshold", "proportion"),
                       C.res, prop, py.nit, en.tol, control) {
    dX <- dim(X)

    if (sum(abs(X[, 1L] - 1)) > .Machine$double.eps^0.75) {
        stop("`X` must have a leading column of 1's")
    }

    ctrl <- initest.control(lambda1 = lambda1,
                            lambda2 = lambda2,
                            numIt = py.nit,
                            eps = en.tol,

                            residCleanMethod = clean.method,
                            residThreshold = C.res,
                            residProportion = prop,
                            pscProportion = prosac,
                            mscaleB = deltaesc,
                            mscaleCC = 1,
                            enpy.control = control)

    ctrl$lambda1 <- ctrl$lambda1 * cc.scale^2
    ctrl$lambda2 <- ctrl$lambda2 * cc.scale^2

    ##
    ## The C++ code needs to now how many observations to *keep*
    ##
    ctrl$pscProportion <- 1 - ctrl$pscProportion

    if (ctrl$residCleanMethod == "proportion") {
        ##
        ## The C++ code needs to now how many observations to *keep*
        ##
        ctrl$residProportion <- 1 - ctrl$residProportion
    }

    ies <- .Call(C_enpy_exact, t(X), y, dX[1L], dX[2L], ctrl)

    return(list(
        coeff = matrix(ies[[1L]], nrow = dX[2L]),
        objF = ies[[2L]]
    ))
}

