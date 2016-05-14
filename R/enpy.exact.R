#' PY (Peña-Yohai) initial estimates for EN
#'
#' Computes the PY initial estimates for EN with exact principal sensitivity components
#'
#' @param X The data matrix X -- a leading column of 1's will be added!
#' @param y The response vector.
#' @param lambda,alpha The EN penalty parameters (NOT adjusted for the number of observations
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
#' @useDynLib penseinit C_enpy_exact C_augtrans
#' @importFrom Rcpp evalCpp
enpy.exact <- function(X, y, alpha, lambda, deltaesc, cc.scale,
                       psc.method = c("Mn", "Qp"), prosac,
                       clean.method = c("threshold", "proportion"),
                       C.res, prop, py.nit, en.tol, control) {
    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X, dX[1L], dX[2L])
    dX[2L] <- dX[2L] + 1L

    ctrl <- initest.control(lambda = lambda,
                            alpha = alpha,
                            numIt = py.nit,
                            eps = en.tol,

                            resid.clean.method = clean.method,
                            resid.threshold = C.res,
                            resid.proportion = prop,
                            psc.proportion = prosac,
                            mscale.delta = deltaesc,
                            mscale.cc = cc.scale,
                            enpy.control = control)

    ##
    ## The C++ code needs to now how many observations to *keep*
    ##
    ctrl$psc.proportion <- 1 - ctrl$psc.proportion

    if (ctrl$resid.clean.method == "proportion") {
        ##
        ## The C++ code needs to now how many observations to *keep*
        ##
        ctrl$resid.proportion <- 1 - ctrl$resid.proportion
    }

    ies <- .Call(C_enpy_exact, Xtr, y, dX[1L], dX[2L], ctrl)

    return(list(
        coeff = matrix(ies[[1L]], nrow = dX[2L]),
        objF = ies[[2L]]
    ))
}

