#' Principal Sensitivity Components
#'
#' Compute the principal sensitiviy components (PSC) for regression.
#'
#' @param X The data matrix X
#' @param y The response vector
#' @param intercept Should an intercept be added or not.
#'
#' @return A numeric matrix with as many rows as \code{X} and as many columns as
#'      PSCs found (at most the number of columns in \code{X} plus one for the intercept).
#'      Each column is a PSC.
#'
#' @references Pena, D., & Yohai, V.. (1999). A Fast Procedure for Outlier Diagnostics in Large
#' Regression Problems. \emph{Journal of the American Statistical Association}, 94(446),
#' 434–445. \url{http://doi.org/10.2307/2670164}
#'
#' @useDynLib penseinit
#' @export
prinsens <- function(X, y, intercept = TRUE) {
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

    ## Add leading column of 1's
    if (identical(intercept, TRUE)) {
        X <- cbind(1, X)
    }

#     numPSCs <- 0L
#     pscres <- .Call("C_pscs", t(X), y, dX[1L], ncol(X), numPSCs,
#                     PACKAGE = "penseinit")
#
#     if (is.null(pscres)) {
#         stop("Could not compute principal sensitivity components. Matrix `x` is singular.")
#     } else if (numPSCs == 0L) {
#         stop("Could not compute principal sensitivity components. All eigenvalues are too small.")
#     }
#
#     pscres <- matrix(pscres[seq_len(dX[1L] * numPSCs)], nrow = dX[1L])

    pscres <- .Call("C_pscs2", t(X), y, dX[1L], ncol(X), PACKAGE = "penseinit")

    if (is.null(pscres)) {
        stop("Could not compute principal sensitivity components. Matrix `x` is singular.")
    } else if (length(pscres) == 0L) {
        stop("Could not compute principal sensitivity components. All eigenvalues are too small.")
    }

    pscres <- matrix(pscres, nrow = dX[1L])

    return(pscres)
}
