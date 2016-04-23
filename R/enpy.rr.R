#' PY (Peña-Yohai) initial estimates for EN
#'
#' Computes the PY initial estimates for EN with approximated principal sensitivity components
#' by the ridge regression solution.
#'
#' @param X The data matrix X -- a leading column of 1's will be added!
#' @param y The response vector.
#' @param lambda1,lambda2 The EN penalty parameters (NOT adjusted for the number of observations
#'          in \code{X}).
#' @param deltaesc,cc.scale Parameters for the M-equation of the scale. The default
#'          rho function is Tukey's bisquare. This can be changed by the parameter \code{control}.
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
#' @useDynLib penseinit C_enpy_rr
#' @importFrom Rcpp evalCpp
enpy.rr <- function(X, y, lambda1, lambda2, deltaesc, cc.scale,
                    prosac, clean.method = c("threshold", "proportion"),
                    C.res, prop, py.nit, en.tol, control) {

    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X, dX[1L], dX[2L])
    dX[2L] <- dX[2L] + 1L

    ctrl <- initest.control(lambda1 = lambda1,
                            lambda2 = lambda2,
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

    usableProp <- 1
    if (ctrl$resid.clean.method == "proportion") {
        ##
        ## The C++ code needs to now how many observations to *keep*
        ##
        ctrl$resid.proportion <- 1 - ctrl$resid.proportion

        usableProp <- ctrl$resid.proportion
    }

    if (ctrl$lambda2 == 0) {
        if (dX[2L] >= dX[1L]) {
            stop("`enpy.rr` can not be used for data with more variables than observations if ",
                 "`lambda2` is 0.")
        } else if (dX[2L] >= ceiling(usableProp * dX[1L])) {
            stop("With the specified proportion of observations to remove, the number of ",
                 "observations will be smaller than the number of variables.\nIn this case ",
                 "`enpy.rr` can not be used with `lambda2` = 0")
        }
    }

    ies <- .Call(C_enpy_rr, Xtr, y, dX[1L], dX[2L], ctrl)

    return(list(
        coeff = matrix(ies[[1L]], nrow = dX[2L]),
        objF = ies[[2L]]
    ))
}

