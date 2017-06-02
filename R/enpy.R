#' PY (Pena-Yohai) initial estimates for the EN S-estimator
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
#' @return \item{initCoef}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
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
    lambda <- .check_arg(lambda, "numeric", range = 0)

    if (alpha < .Machine$double.eps) {
        options$pscMethod <- 'rr'
    }

    result <- switch(
        options$pscMethod,
        rr = enpy_rr(X, y, alpha, lambda, options, en_options),
        exact = enpy_exact(X, y, alpha, lambda, options, en_options)
    )

    resorder <- sort.list(result$objF, na.last = NA, method = "quick")

    dups <- which(diff(result$objF[resorder]) < en.tol)
    if (length(dups) > 0) {
        remove <- resorder[dups]
        result$objF <- result$objF[-remove]
        result$coeff <- result$coeff[ , -remove, drop = FALSE]
    }

    return(result)
}
