#' PY (Peña-Yohai) initial estimates for S-estimates of regression
#'
#' Computes the PY initial estimates for S-estimates of regression
#'
#' @param X The data matrix X
#' @param y The response vector
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
#'
#' @return \item{initCoef}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
#'
#' @useDynLib penseinit
#' @export
pyinit <- function(X, y, deltasc, cc.scale, prosac,
                 clean.method = c("threshold", "proportion"), C.res, prop,
                 py.nit, en.tol) {
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

    ctrl <- initest.control(numIt = py.nit,
                            eps = en.tol,
                            residCleanMethod = clean.method,
                            residThreshold = C.res,
                            residProportion = prop,
                            pscProportion = prosac,

                            mscaleB = deltasc,
                            mscaleCC = 1,

                            ## We don't need those parameters
                            lambda1 = 0,
                            lambda2 = 0)

    cnuminits <- 0L

    ies <- .Call("C_initpy", t(X), y, nrow(X), ncol(X),
                 ctrl$numIt,
                 ctrl$eps,
                 ctrl$residThreshold,
                 ctrl$residProportion,
                 ctrl$pscProportion,
                 ctrl$mscaleB,
                 ctrl$mscaleCC,
                 ctrl$mscaleMaxIt,
                 ctrl$mscaleEPS,
                 ctrl$mscaleRhoFun,
                 cnuminits, PACKAGE = "penseinit")

    ies <- matrix(ies[cnuminits * (dX[2L] + 1L)], nrow = dX[2L] + 1L, ncol = cnuminits)

    return(list(
        initCoef = ies,
        objF = numeric(cnuminits)
    ))
}
