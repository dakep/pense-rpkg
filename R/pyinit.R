#' PY (Pena-Yohai) initial estimates for S-estimates of regression
#'
#' Computes the PY initial estimates for S-estimates of regression.
#'
#' @param X the data matrix X.
#' @param y the response vector.
#' @param intercept should an intercept be included in the models. Defaults to \code{TRUE}.
#' @param deltaesc,cc.scale parameters for the M-equation of the scale. If \code{cc.scale}
#'          is missing or invalid, it will be chosen such that the expected value of the
#'          rho function under the normal model is equal to \code{deltaesc}. To specify the
#'          rho function, see parameter \code{control}.
#' @param prosac the proportion of observations to remove based on PSCs.
#' @param clean.method how to clean the data based on large residuals.
#'          If \code{"threshold"}, all observations with scaled residuals larger than
#'          \code{C.res} will be removed (\code{C.res} corresponds to the constant
#'          \eqn{C_1} from equation (21) in Pena & Yohai (1999).
#'          If \code{"proportion"}, observations with the largest \code{prop} residuals
#'          will be removed.
#' @param py.nit the maximum number of iterations to perform.
#' @param en.tol the relative tolerance for convergence.
#' @param control optional further control parameters from \code{\link{enpy.control}}.
#' @param C.res,prop see parameter \code{clean.method} for details.
#'
#' @return \item{initCoef}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
#'
#' @references Pena, D., & Yohai, V.. (1999). A Fast Procedure for Outlier Diagnostics in Large
#' Regression Problems. \emph{Journal of the American Statistical Association}, 94(446),
#' 434-445. \url{http://doi.org/10.2307/2670164}
#'
#' @useDynLib pense C_initpy
#' @export
pyinit <- function(X, y, intercept = TRUE, deltaesc, cc.scale, prosac,
                   clean.method = c("threshold", "proportion"), C.res, prop,
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

    if (!identical(intercept, FALSE)) {
        ## Add leading column of 1's
        X <- cbind(1, X)
    }

    ctrl <- initest.control(numIt = py.nit,
                            eps = en.tol,
                            resid.clean.method = clean.method,
                            resid.threshold = C.res,
                            resid.proportion = prop,
                            psc.proportion = prosac,

                            mscale.delta = deltaesc,
                            mscale.cc = cc.scale,
                            enpy.control = control,

                            ## We don't need those parameters
                            lambda = 0,
                            alpha = 0)

    ##
    ## The C code needs to now how many observations to *keep*
    ##
    ctrl$psc.proportion <- 1 - ctrl$psc.proportion

    usableProp <- 1
    if (ctrl$resid.clean.method == "proportion") {
        ## The C code needs to now how many observations to *keep*
        ctrl$resid.proportion <- 1 - ctrl$resid.proportion

        usableProp <- ctrl$resid.proportion
    }

    if (dX[2L] >= dX[1L]) {
        stop("`pyinit` can not be used for data with more variables than observations")
    } else if (dX[2L] >= ceiling(usableProp * dX[1L])) {
        stop("With the specified proportion of observations to remove, the number of ",
             "observations will be smaller than the number of variables.\nIn this case ",
             "`pyinit` can not be used.")
    }

    dX <- dim(X)

    ies <- .Call(C_initpy, t(X), y, dX[1L], dX[2L],
                 ctrl$numIt,
                 ctrl$eps,
                 ctrl$resid.threshold,
                 ctrl$resid.proportion,
                 ctrl$psc.proportion,
                 ctrl$mscale.delta,
                 ctrl$mscale.cc,
                 ctrl$mscale.maxit,
                 ctrl$mscale.tol,
                 ctrl$mscale.rho.fun)

    initCoefs <- matrix(ies[[2L]][seq_len(ies[[1L]] * dX[2L])], nrow = dX[2L], ncol = ies[[1L]])

    return(list(
        initCoef = initCoefs,
        objF = ies[[3L]][seq_len(ies[[1L]])]
    ))
}
