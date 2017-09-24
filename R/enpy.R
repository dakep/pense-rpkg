#' PY (Pena-Yohai) initial estimates for EN S-estimators
#'
#' Computes the PY initial estimates for the EN S-estimator with different
#' strategies for computing the principal sensitivity components.
#'
#' Two different methods to calculate the sensitivity components are implemented:
#' \describe{
#'      \item{\code{"rr"}}{Approximate the PSCs by using the residuals from the
#'          elastic net fit and the hat matrix from the ridge regression.
#'          This method only works if \code{alpha} < 1 or
#'          \code{ncol(X)} < \code{nrow(X)}.}
#'      \item{\code{"exact"}}{Calculate the PSCs from the difference between the
#'          residuals and leave-one-out residuals from elastic net.}
#' }
#'
#' @param X data matrix with predictors.
#' @param y response vector.
#' @param alpha,lambda EN penalty parameters (NOT adjusted for the number of
#'      observations in \code{X}).
#' @param options additional options for the initial estimator. See
#'      \code{\link{initest_options}} for details.
#' @param en_options additional options for the EN algorithm. See
#'      \code{\link{en_options}} for details.
#'
#' @return \item{coeff}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
#'
#' @references Pena, D., and Yohai, V.J. (1999).
#'     A Fast Procedure for Outlier Diagnostics in Large Regression Problems.
#'     \emph{Journal of the American Statistical Association}, \bold{94}(446),
#'     434-445. \url{http://doi.org/10.2307/2670164}
#'
#' @example examples/enpy.R
#'
#' @export
enpy <- function(X, y, alpha, lambda,
                 options = initest_options(),
                 en_options = en_options_aug_lars()) {
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

    alpha <- .check_arg(alpha, "numeric", range = c(0, 1),
                        range_test_lower = ">=", range_test_upper = "<=")
    lambda <- .check_arg(lambda, "numeric", range = 0, range_test_lower = ">=")

    if (alpha < .Machine$double.eps) {
        options$pscMethod <- 'rr'
    }

    if (lambda == 0) {
        options$pscMethod <- "ols"
    }

    result <- switch(
        options$pscMethod,
        ols = enpy_ols(X, y, options),
        rr = enpy_rr(X, y, alpha, lambda, options, en_options),
        exact = enpy_exact(X, y, alpha, lambda, options, en_options)
    )

    resorder <- sort.list(result$objF, na.last = NA, method = "quick")

    dups <- which(diff(result$objF[resorder]) < options$eps)
    if (length(dups) > 0) {
        remove <- resorder[dups]
        result$objF <- result$objF[-remove]
        result$coeff <- result$coeff[ , -remove, drop = FALSE]
    }

    return(result)
}


## PY (Pena-Yohai) initial estimates for EN S-Estimators
##
## Computes the PY initial estimates for EN  S-Estimators with exact
## principal sensitivity components.
##
## @param X data matrix with predictors -- a leading column of 1's will be
##      added!
## @param y response vector.
## @param lambda,alpha The EN penalty parameters (NOT adjusted for the number
##      of observations in \code{X}).
## @param options additional options for the initial estimator. See
##      \code{\link{initest_options}} for details.
## @param en_options additional options for the EN algorithm. See
##      \code{\link{en_options}} for details.
##
## @return \item{coeff}{A numeric matrix with one initial coefficient per
##      column}
##      \item{objF}{A vector of values of the objective function for the
##      respective coefficient}
##
#' @useDynLib pense, .registration = TRUE
#' @importFrom Rcpp evalCpp
enpy_exact <- function(X, y, alpha, lambda, options, en_options) {
    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X)
    dX[2L] <- dX[2L] + 1L

    ies <- .Call(C_enpy_exact, Xtr, y, alpha, lambda, options, en_options)

    return(list(
        coeff = matrix(ies[[1L]], nrow = dX[2L]),
        objF = ies[[2L]]
    ))
}

## PY (Pena-Yohai) initial estimates for EN S-Estimators
##
## Computes the PY initial estimates for EN S-Estimators with
## principal sensitivity components approximated by the ridge regression
## solution.
##
## @param X data matrix with predictors -- a leading column of 1's will be
##      added!
## @param y response vector.
## @param lambda,alpha The EN penalty parameters (NOT adjusted for the number
##      of observations in \code{X}).
## @param options additional options for the initial estimator. See
##      \code{\link{initest_options}} for details.
## @param en_options additional options for the EN algorithm. See
##      \code{\link{en_options}} for details.
##
## @return \item{coeff}{A numeric matrix with one initial coefficient per column}
##         \item{objF}{A vector of values of the objective function for the respective coefficient}
##
#' @useDynLib pense, .registration = TRUE
#' @importFrom Rcpp evalCpp
enpy_rr <- function(X, y, alpha, lambda, options, en_options) {
    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X)
    dX[2L] <- dX[2L] + 1L

    usableProp <- with(options, if(keepResidualsMethod == 1) {
        keepResidualsProportion * keepPSCProportion
    } else {
        keepPSCProportion
    })

    if (alpha >= 1 - .Machine$double.eps) {
        if (dX[2L] >= dX[1L]) {
            stop("`enpy_rr` can not be used for data with more variables than ",
                 "observations if `alpha` is 1.")
        } else if (dX[2L] >= ceiling(usableProp * dX[1L])) {
            stop("With the specified proportion of observations to remove, ",
                 "the number of observations will be smaller than the number ",
                 "of variables.\n",
                 "In this case `enpy_rr` can not be used when `alpha` is 1")
        }
    }

    ies <- .Call(C_enpy_rr, Xtr, y, alpha, lambda, options, en_options)

    return(list(
        coeff = matrix(ies[[1L]], nrow = dX[2L]),
        objF = ies[[2L]]
    ))
}


## PY (Pena-Yohai) initial estimates for EN S-Estimators
##
## Computes the PY initial estimates for EN S-Estimators with
## principal sensitivity components approximated by the ridge regression
## solution.
##
## @param X data matrix with predictors -- a leading column of 1's will be
##      added!
## @param y response vector.
## @param options additional options for the initial estimator. See
##      \code{\link{initest_options}} for details.
##
## @return \item{coeff}{A numeric matrix with one initial coefficient per column}
##         \item{objF}{A vector of values of the objective function for the respective coefficient}
##
#' @useDynLib pense, .registration = TRUE
#' @importFrom Rcpp evalCpp
enpy_ols <- function(X, y, options) {
    dX <- dim(X)

    Xtr <- .Call(C_augtrans, X)
    dX[2L] <- dX[2L] + 1L

    ies <- .Call(C_py_ols, Xtr, y, options)

    return(list(
        coeff = matrix(ies[[1L]], nrow = dX[2L]),
        objF = ies[[2L]]
    ))
}
