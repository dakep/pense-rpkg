#' Principal Sensitivity Components
#'
#' Compute the principal sensitiviy components (PSC) for regression.
#'
#' @param X The data matrix X
#' @param y The response vector
#' @param method Use ordinary least squares (\code{"ols"}) or elastic net (\code{"en"}) to compute
#'               the PSCs.
#' @param intercept Should an intercept be added or not.
#' @param alpha,lambda The values for the parameters controlling the penalization for elastic net.
#' @param maxit The maximum number of iterations for elastic net.
#' @param eps The relative tolerance for convergence for gradient-descent (default 1e-8) or
#'      the threshold for treating numbers as 0 in LARS (default .Machine$double.eps)
#' @param centering Should the rows be centered first for elastic net.
#' @param en.algorithm algorithm to use to compute the elastic net solution.
#'
#' @return A numeric matrix with as many rows as \code{X} and as many columns as
#'      PSCs found (at most the number of columns in \code{X} plus one for the intercept).
#'      Each column is a PSC.
#'
#' @references Pena, D., & Yohai, V.. (1999). A Fast Procedure for Outlier Diagnostics in Large
#' Regression Problems. \emph{Journal of the American Statistical Association}, 94(446),
#' 434-445. \url{http://doi.org/10.2307/2670164}
#'
#' @useDynLib pense C_pscs_ols C_pscs_en C_augtrans
#' @export
prinsens <- function(X, y, method = c("ols", "en"), intercept = TRUE,
                     alpha, lambda, maxit = 50000, eps, centering = TRUE,
                     en.algorithm = c("augmented-lars",
                                      "coordinate-descent",
                                      "augmented-lars-gram",
                                      "augmented-lars-nogram")) {
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
        Xtr <- .Call(C_augtrans, X, dX[1L], dX[2L])
        dX[2L] <- dX[2L] + 1L
    } else {
        Xtr <- t(X)
    }

    method <- match.arg(method)

    if (method == "en") {
        if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) || alpha < 0 || alpha > 1) {
            stop("`alpha` must be single number between in the range [0, 1]")
        }

        if (length(lambda) != 1L || !is.numeric(lambda) || is.na(lambda) || lambda < 0) {
            stop("`lambda` must be single number >= 0")
        }

        if (length(maxit) != 1L || !is.numeric(maxit) || is.na(maxit) || maxit <= 1) {
            stop("`maxit` must be single integer > 1")
        }

        en.algorithm <- match.arg(en.algorithm);
        if (missing(eps)) {
            eps <- switch(en.algorithm, `coordinate-descent` = 1e-8, .Machine$double.eps)
        }

        if (length(eps) != 1L || !is.numeric(eps) || is.na(eps) || eps <= 0) {
            stop("`eps` must be single number > 0")
        }

        if (length(centering) != 1L || !is.logical(centering) || is.na(centering)) {
            warning("`centering` must be single logical value. Using TRUE as default.")
        }

        alpha <- as.numeric(alpha)
        lambda <- as.numeric(lambda)
        maxit <- as.integer(maxit)
        centering <- 1L - as.integer(identical(centering, FALSE))
        en.algorithm <- .enalgo2IntEnalgo(en.algorithm)

        pscres <- .Call(C_pscs_en, Xtr, y, dX[1L], dX[2L],
                        alpha, lambda, maxit, eps, centering, en.algorithm)

        if (is.null(pscres)) {
            stop("Could not compute principal sensitivity components.")
        }
    } else {
#         numPSCs <- 0L
#         pscres <- .Call("C_pscs", t(X), y, dX[1L], ncol(X), numPSCs,
#                         PACKAGE = "pense")
#
#         if (is.null(pscres)) {
#             stop("Could not compute principal sensitivity components. Matrix `x` is singular.")
#         } else if (numPSCs == 0L) {
#             stop("Could not compute principal sensitivity components. All eigenvalues are too small.")
#         }
#
#         pscres <- matrix(pscres[seq_len(dX[1L] * numPSCs)], nrow = dX[1L])

        pscres <- .Call(C_pscs_ols, Xtr, y, dX[1L], dX[2L])

        if (is.null(pscres)) {
            stop("Could not compute principal sensitivity components. Matrix `x` is singular.")
        }
    }

    if (length(pscres) == 0L) {
        stop("Could not compute principal sensitivity components. All eigenvalues are too small.")
    }

    pscres <- matrix(pscres, nrow = dX[1L])

    return(pscres)
}
