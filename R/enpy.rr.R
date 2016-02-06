#' PY (Peña-Yohai) initial estimates for EN
#'
#' Computes the PY initial estimates for EN with approximated principal sensitivity components
#' by the ridge regression solution.
#'
#' @param X The data matrix X, with leading column of 1's
#' @param y The response vector
#' @param lambda1,lambda2 The EN penalty parameters (adjusted for the number of observations
#'          in \code{X})
#' @param deltasc,cc.scale Parameters for the M-equation of the scale. Tukey's bisquare rho function
#'      is used internally.
#' @param prosac The proportion of observations to remove based on PSCs
#' @param clean.method How to clean the data based on large residuals.
#'          If \code{"threshold"}, all observations with scaled residuals larger than
#'          \code{C.res} will be removed.
#'          If \code{"proportion"}, observations with the largest \code{prop} residuals
#'          will be removed.
#' @param py.nit The maximum number of iterations to perform.
#' @param en.tol The relative tolerance for convergence.
#' @param en.centering Should rows with a leading 1 be centered in the elastic net algorithm.
#'      Default's to \code{TRUE}.
#'
#' @return \item{initCoef}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
#'
#' @useDynLib penseinit
#' @importFrom Rcpp evalCpp
enpy.rr <- function(X, y, lambda1, lambda2, deltasc, cc.scale,
                    prosac, clean.method = c("threshold", "proportion"),
                    C.res, prop, py.nit, en.tol, en.centering = TRUE) {

    dX <- dim(X)

    if (sum(abs(X[, 1L] - 1)) > 1e-15) {
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

                            en.centering = en.centering,

                            mscaleB = deltasc,
                            mscaleCC = 1)

    ctrl$lambda1 <- ctrl$lambda1 * cc.scale^2

    usableProp <- ctrl$pscProportion
    if (ctrl$residCleanMethod == "proportion") {
        usableProp <- ctrl$residProportion * ctrl$pscProportion
    }

    if (ctrl$lambda2 == 0) {
        if (dX[2L] >= dX[1L]) {
            stop("`enpy.rr` can not be used for data with more variables than observations if ",
                 "`lambda2` is 0.")
        } else if (dX[2L] >= as.integer(usableProp * dX[1L])) {
            stop("With the specified proportion of observations to remove, the number of ",
                 "observations will be smaller than the number of variables.\nIn this case ",
                 "`enpy.rr` can not be used with `lambda2` = 0")
        }
    }

    ies <- .Call("C_enpy_rr", t(X), y, dX[1L], dX[2L], ctrl, PACKAGE = "penseinit")

    return(list(
        initCoef = matrix(ies[[1L]], nrow = dX[2L]),
        objF = ies[[2L]]
    ))
}

