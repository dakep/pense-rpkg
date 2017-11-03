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
#'          \code{ncol(x)} < \code{nrow(x)}.}
#'      \item{\code{"exact"}}{Calculate the PSCs from the difference between the
#'          residuals and leave-one-out residuals from elastic net.}
#' }
#'
#' @param x data matrix with predictors.
#' @param y response vector.
#' @param alpha,lambda EN penalty parameters (NOT adjusted for the number of
#'      observations in \code{x}).
#' @param delta desired breakdown point of the resulting estimator.
#' @param cc tuning constant for the S-estimator. Default is to chosen based
#'      on the breakdown point \code{delta}. Should never have to be changed.
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
enpy <- function(x, y, alpha, lambda, delta, cc,
                 options = initest_options(),
                 en_options = en_options_aug_lars()) {
    y <- drop(y)

    dx <- dim(x)
    dY <- dim(y)
    yl <- length(y)

    if (is.null(yl) || (!is.null(dY) && length(dY) != 1L) || !is.numeric(y)) {
        stop("`yl` must be a numeric vector")
    }

    if (is.null(dx) || length(dx) != 2L || !is.numeric(x) || dx[1L] != yl) {
        stop("`x` must be a numeric matrix with the same number of observations as `y`")
    }

    if (anyNA(x) || anyNA(y)) {
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

    options$mscaleDelta <- .check_arg(
        delta,
        "numeric",
        range = c(0, 0.5),
        range_test_upper = "<="
    )

    options$mscaleCC <- if (missing(cc)) {
        consistency.rho(options$mscaleDelta, 1L)
    } else {
        .check_arg(
            cc,
            "numeric",
            range = 0
        )
    }

    result <- switch(
        options$pscMethod,
        ols = enpy_ols(x, y, options),
        rr = enpy_rr(x, y, alpha, lambda, options, en_options),
        exact = enpy_exact(x, y, alpha, lambda, options, en_options)
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
## @param x data matrix with predictors -- a leading column of 1's will be
##      added!
## @param y response vector.
## @param lambda,alpha The EN penalty parameters (NOT adjusted for the number
##      of observations in \code{x}).
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
enpy_exact <- function(x, y, alpha, lambda, options, en_options) {
    dx <- dim(x)

    xtr <- .Call(C_augtrans, x)
    dx[2L] <- dx[2L] + 1L

    ies <- .Call(C_enpy_exact, xtr, y, alpha, lambda, options, en_options)

    return(list(
        coeff = matrix(ies[[1L]], nrow = dx[2L]),
        objF = ies[[2L]]
    ))
}

## PY (Pena-Yohai) initial estimates for EN S-Estimators
##
## Computes the PY initial estimates for EN S-Estimators with
## principal sensitivity components approximated by the ridge regression
## solution.
##
## @param x data matrix with predictors -- a leading column of 1's will be
##      added!
## @param y response vector.
## @param lambda,alpha The EN penalty parameters (NOT adjusted for the number
##      of observations in \code{x}).
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
enpy_rr <- function(x, y, alpha, lambda, options, en_options) {
    dx <- dim(x)

    xtr <- .Call(C_augtrans, x)
    dx[2L] <- dx[2L] + 1L

    usableProp <- with(options, if(keepResidualsMethod == 1) {
        keepResidualsProportion * keepPSCProportion
    } else {
        keepPSCProportion
    })

    if (alpha >= 1 - .Machine$double.eps) {
        if (dx[2L] >= dx[1L]) {
            stop("`enpy_rr` can not be used for data with more variables than ",
                 "observations if `alpha` is 1.")
        } else if (dx[2L] >= ceiling(usableProp * dx[1L])) {
            stop("With the specified proportion of observations to remove, ",
                 "the number of observations will be smaller than the number ",
                 "of variables.\n",
                 "In this case `enpy_rr` can not be used when `alpha` is 1")
        }
    }

    ies <- .Call(C_enpy_rr, xtr, y, alpha, lambda, options, en_options)

    return(list(
        coeff = matrix(ies[[1L]], nrow = dx[2L]),
        objF = ies[[2L]]
    ))
}


## PY (Pena-Yohai) initial estimates for EN S-Estimators
##
## Computes the PY initial estimates for EN S-Estimators with
## principal sensitivity components approximated by the ridge regression
## solution.
##
## @param x data matrix with predictors -- a leading column of 1's will be
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
enpy_ols <- function(x, y, options) {
    dx <- dim(x)

    xtr <- .Call(C_augtrans, x)
    dx[2L] <- dx[2L] + 1L

    ies <- .Call(C_py_ols, xtr, y, options)

    return(list(
        coeff = matrix(ies[[1L]], nrow = dx[2L]),
        objF = ies[[2L]]
    ))
}
