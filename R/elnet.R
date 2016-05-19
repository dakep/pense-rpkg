#' Elastic Net Regression
#'
#' Compute the elastic net regression coefficients
#'
#' This solves the minimization problem
#' \deqn{\frac{1}{2 N} RSS + \lambda \left( \frac{(1 - \alpha)} {2} \| \beta \|_2^2 + \alpha \| \beta \|_1 \right)}{
#'      (1/2N) RSS + \lambda * ( (1 - \alpha) / 2 * L2(\beta)^2 + \alpha * L1(\beta) )}
#'
#' If weights are supplied, the minimization problem becomes
#' \deqn{\frac{1}{2 N} \sum_{i = 1}^n w_i r_i^2 + \lambda \left( \frac{(1 - \alpha)} {2} \| \beta \|_2^2 + \alpha \| \beta \|_1 \right)}{
#'      (1/2N) sum(w * (y - r^2) + \lambda * ( (1 - \alpha) / 2 * L2(\beta)^2 + \alpha * L1(\beta) )}
#'
#'
#' @param X The data matrix X
#' @param y The response vector
#' @param alpha,lambda The values for the parameters controlling the penalization
#' @param weights an optional vector of weights to be used in the fitting process. Should be
#'                \code{NULL} or a numeric vector. If non-NULL, weighted EN is used with weights
#'                \code{weights} See also ‘Details’.
#' @param maxit The maximum number of iterations
#' @param eps The relative tolerance for convergence for gradient-descent (default 1e-8) or
#'      the threshold for treating numbers as 0 in LARS (default .Machine$double.eps)
#' @param centering Should the rows be centered first
#' @param addLeading1s Should a leading column of 1's be appended? If \code{FALSE}, this has
#'      to be done before calling this function.
#' @param en.algorithm algorithm to use to compute the elastic net solution.
#'
#' @return \item{coefficients}{The regression coefficients}
#'         \item{residuals}{The residuals}
#'         \item{converged}{Did the algorithm converge?}
#'
#' @useDynLib penseinit C_elnet
#' @export
elnet <- function(X, y, alpha, lambda, weights, maxit = 10000, eps, centering = TRUE,
                  addLeading1s = TRUE, en.algorithm = c("augmented-lars",
                                                        "coordinate-descent",
                                                        "augmented-lars-gram",
                                                        "augmented-lars-nogram")) {
    y <- drop(y)

    dX <- dim(X)
    dY <- dim(y)
    yl <- length(y)

    weighted <- FALSE

    if (is.null(yl) || (!is.null(dY) && length(dY) != 1L) || !is.numeric(y)) {
        stop("`yl` must be a numeric vector")
    }

    if (is.null(dX) || length(dX) != 2L || !is.numeric(X) || dX[1L] != yl) {
        stop("`X` must be a numeric matrix with the same number of observations as `y`")
    }

    if (missing(weights)) {
        weights <- NULL
    }

    if (!is.null(weights)) {
        weighted <- TRUE
        weights <- drop(weights)

        if (length(weights) != dX[1L] || !is.null(dim(weights))) {
            stop("If `weights` are supplied, they must be a numeric vector the same length as `y`")
        }
    }

    if (anyNA(X) || anyNA(y) || anyNA(weights)) {
        stop("Missing values are not supported")
    }


    if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) || alpha < 0 || alpha > 1) {
        stop("`alpha` must be single number in the range [0, 1]")
    }

    if (length(lambda) != 1L || !is.numeric(lambda) || is.na(lambda) || lambda < 0) {
        stop("`lambda` must be single number >= 0")
    }

    if (length(maxit) != 1L || !is.numeric(maxit) || is.na(maxit) || maxit <= 1) {
        stop("`maxit` must be single integer > 1")
    }

    en.algorithm <- match.arg(en.algorithm)

    if (missing(eps)) {
        eps <- switch(en.algorithm, `coordinate-descent` = 1e-8, .Machine$double.eps)
    }

    if (length(eps) != 1L || !is.numeric(eps) || is.na(eps) || eps <= 0) {
        stop("`eps` must be single number > 0")
    }

    if (length(centering) != 1L || !is.logical(centering) || is.na(centering)) {
        warning("`centering` must be single logical value. Using TRUE as default.")
    }

    if (weighted) {
        elnetres <- .elnet.wfit(X, y, weights, alpha, lambda, maxit, eps, centering, addLeading1s,
                                warmCoef = NULL, en.algorithm = en.algorithm)
    } else {
        elnetres <- .elnet.fit(X, y, alpha, lambda, maxit, eps, centering, addLeading1s,
                               warmCoef = NULL, en.algorithm = en.algorithm)
    }

    return(elnetres)
}

## Internal function to fit an EN linear regression WITHOUT parameter checks!
#' @useDynLib penseinit C_augtrans C_elnet
.elnet.fit <- function(X, y, alpha, lambda, maxit, eps, centering = TRUE, addLeading1s = TRUE,
                       warmCoef = NULL, en.algorithm) {
    y <- drop(y)
    dX <- dim(X)

    ## Add leading column of 1's
    if (!identical(addLeading1s, FALSE)) {
        Xtr <- .Call(C_augtrans, X, dX[1L], dX[2L])
    } else {
        Xtr <- t(X)
    }

    warm <- 0L

    if (!is.null(warmCoef)) {
        warm <- 1L
    }

    alpha <- as.numeric(alpha)
    lambda <- as.numeric(lambda)
    maxit <- as.integer(maxit)
    centering <- 1L - as.integer(identical(centering, FALSE))
    en.algorithm <- .enalgo2IntEnalgo(en.algorithm)

    elnetres <- .Call(C_elnet, Xtr, y,
                      warmCoef,
                      dX[1L], nrow(Xtr),
                      alpha,
                      lambda,
                      maxit,
                      eps,
                      centering,
                      warm,
                      en.algorithm)

    if (!identical(elnetres[[1L]], TRUE)) {
        warning("Elastic Net algorithm did not converge.")
    }

    names(elnetres) <- c("converged", "coefficients", "residuals")

    return(elnetres)
}

## Internal function to fit an EN linear regression WITHOUT parameter checks!
#' @useDynLib penseinit C_augtrans C_elnet_weighted
.elnet.wfit <- function(X, y, weights, alpha, lambda, maxit, eps, centering = TRUE,
                        addLeading1s = TRUE, warmCoef = NULL, en.algorithm) {
    y <- drop(y)
    dX <- dim(X)

    ## Add leading column of 1's
    if (!identical(addLeading1s, FALSE)) {
        Xtr <- .Call(C_augtrans, X, dX[1L], dX[2L])
    } else {
        Xtr <- t(X)
    }

    warm <- 0L

    if (!is.null(warmCoef)) {
        warm <- 1L
    }

    weights <- as.numeric(weights)
    alpha <- as.numeric(alpha)
    lambda <- as.numeric(lambda)
    maxit <- as.integer(maxit)
    centering <- 1L - as.integer(identical(centering, FALSE))
    en.algorithm <- .enalgo2IntEnalgo(en.algorithm)

    elnetres <- .Call(C_elnet_weighted, Xtr, y,
                      weights,
                      warmCoef,
                      dX[1L], nrow(Xtr),
                      alpha,
                      lambda,
                      maxit,
                      eps,
                      centering,
                      warm,
                      en.algorithm)

    if (!identical(elnetres[[1L]], TRUE)) {
        warning("Weighted Elastic Net algorithm did not converge.")
    }

    names(elnetres) <- c("converged", "coefficients", "residuals")

    return(elnetres)
}
