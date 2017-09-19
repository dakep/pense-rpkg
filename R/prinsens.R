#' Principal Sensitivity Components
#'
#' Compute the principal sensitiviy components (PSC) for regression.
#'
#' @param X data matrix with predictors
#' @param y response vector
#' @param method use ordinary least squares (\code{"ols"}) or elastic net
#'      (\code{"en"}) to compute the PSCs.
#' @param intercept Should an intercept be added or not.
#' @param alpha,lambda The values for the parameters controlling the
#'      penalization for elastic net.
#' @param en_options additional options for the EN algorithm. See
#'      \code{\link{en_options}} for details.
#'
#' @return A numeric matrix with as many rows as \code{X} and as many columns as
#'      PSCs found (at most the number of columns in \code{X} plus one for the
#'      intercept). Each column is a PSC.
#'
#' @references Pena, D., and Yohai, V.J. (1999).
#'     A Fast Procedure for Outlier Diagnostics in Large Regression Problems.
#'     \emph{Journal of the American Statistical Association}, \bold{94}(446),
#'     434-445. \url{http://doi.org/10.2307/2670164}
#'
#' @useDynLib pense, .registration = TRUE
#' @export
prinsens <- function(X, y, method = c("ols", "en"), intercept = TRUE,
                     alpha, lambda, en_options = en_options_aug_lars()) {
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

    intercept <- .check_arg(intercept, "logical")

    ## Add leading column of 1's
    Xtr <- if (isTRUE(intercept)) {
        .Call(C_augtrans, X)
    } else {
        t(X)
    }

    method <- match.arg(method)

    if (method == "en") {
        alpha <- .check_arg(alpha, "numeric", range = c(0, 1),
                            range_test_lower = ">=", range_test_upper = "<=")
        lambda <- .check_arg(lambda, "numeric", range = 0)

        pscres <- .Call(
            C_pscs_en,
            Xtr,
            y,
            alpha,
            lambda,
            intercept,
            en_options
        )

        if (is.null(pscres)) {
            stop("Could not compute principal sensitivity components.")
        }
    } else {
        pscres <- .Call(C_pscs_ols, Xtr, y)

        if (is.null(pscres)) {
            stop("Could not compute principal sensitivity components. ",
                 "Matrix `x` is singular.")
        }
    }

    if (length(pscres) == 0L) {
        stop("Could not compute principal sensitivity components.",
             "All eigenvalues are too small.")
    }

    pscres <- matrix(pscres, nrow = dX[1L])

    return(pscres)
}
