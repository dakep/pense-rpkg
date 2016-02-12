#' Elastic Net Regression
#'
#' Compute the elastic net regression coefficients
#'
#' This solves the minimization problem
#' \deqn{\frac{1}{2 N} RSS + \lambda \left( (1 - \alpha) \frac{1}{2} \| \beta \|_1 + \alpha \| \beta \|_2^2 \right)}{
#'      (1/2N) RSS + \lambda ( (1 - \alpha) (1/2) L1(\beta) + \alpha L2(\beta)^2 )}
#'
#' @param X The data matrix X
#' @param y The response vector
#' @param alpha,lambda The values for the parameters controlling the penalization
#' @param maxit The maximum number of iterations
#' @param eps The relative tolerance for convergence
#' @param centering Should the rows be centered first
#' @param addLeading1s Should a leading column of 1's be appended? If \code{FALSE}, this has
#'      to be done before calling this function.
#'
#' @return \item{coefficients}{The regression coefficients}
#'         \item{residuals}{The residuals}
#'         \item{converged}{Did the algorithm converge?}
#'
#' @useDynLib penseinit
#' @export
elnet <- function(X, y, alpha, lambda, maxit = 500, eps = 1e-6, centering = TRUE,
                  addLeading1s = TRUE) {
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

    if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) || alpha < 0 || alpha > 1) {
        stop("`alpha` must be single number between in the range [0, 1]")
    }

    if (length(lambda) != 1L || !is.numeric(lambda) || is.na(lambda) || lambda < 0) {
        stop("`lambda` must be single number >= 0")
    }

    if (length(maxit) != 1L || !is.numeric(maxit) || is.na(maxit) || maxit <= 1) {
        stop("`maxit` must be single integer > 1")
    }

    if (length(eps) != 1L || !is.numeric(eps) || is.na(eps) || eps <= 0) {
        stop("`eps` must be single number > 0")
    }

    if (length(centering) != 1L || !is.logical(centering) || is.na(centering)) {
        warning("`centering` must be single logical value. Using TRUE as default.")
    }


    ## Add leading column of 1's
    if (!identical(addLeading1s, FALSE)) {
        X <- cbind(1, X)
    }

    alpha <- as.numeric(alpha)
    lambda <- as.numeric(lambda)
    maxit <- as.integer(maxit)
    centering <- 1L - as.integer(identical(centering, FALSE))

    elnetres <- .Call("C_elnet", t(X), y, dX[1L], ncol(X),
                      alpha,
                      lambda,
                      maxit,
                      eps,
                      centering,
                      PACKAGE = "penseinit")

    if (!identical(elnetres[[1L]], TRUE)) {
        warning("Elastic Net algorithm did not converge.")
    }

    names(elnetres) <- c("converged", "coefficients", "residuals")

    return(elnetres)
}
