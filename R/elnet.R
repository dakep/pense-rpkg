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
#' @param X data matrix with predictors
#' @param y response vector
#' @param alpha,lambda values for the parameters controlling the penalization
#' @param weights an optional vector of weights to be used in the fitting
#'      process. Should be \code{NULL} or a numeric vector. If non-NULL,
#'      weighted EN is used with weights \code{weights}. See also 'Details'.
#' @param options additional options for the EN algorithm. See
#'      \code{\link{en_options}} for details.
#' @param intercept should an intercept be estimated?
#' @param addLeading1s should a leading column of 1's be appended?
#'      If \code{FALSE}, this has to be done before calling this function.
#'
#' @return \item{coefficients}{The regression coefficients}
#'         \item{residuals}{The residuals}
#'         \item{converged}{Did the algorithm converge?}
#'
#' @export
elnet <- function(X, y, alpha, lambda, weights, intercept = TRUE,
                  addLeading1s = TRUE, options = en_options_aug_lars()) {
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
        weights <- .check_arg(weights, "numeric", range = 0, length = yl,
                              range_test_lower = ">=")
    }

    if (anyNA(X) || anyNA(y) || anyNA(weights)) {
        stop("Missing values are not supported")
    }

    intercept <- isTRUE(intercept)

    alpha <- .check_arg(alpha, "numeric", range = c(0, 1),
                        range_test_lower = ">=", range_test_upper = "<=")
    lambda <- .check_arg(lambda, "numeric", range = 0, length = NULL)
    lambda <- sort(lambda, decreasing = TRUE)

    # Check the size of X and y
    if (identical(dX[1L], 0L)) {
        # no observations
        return(list(
            status = 0L,
            message = "no observations",
            coefficients = matrix(NA_real_, nrow = dX[2L] + intercept,
                                  ncol = length(lambda)),
            residuals = matrix(NA_real_, nrow = 0L, ncol = length(lambda))
        ))
    } else if (identical(dX[2L], 0L)) {
        # no predictors given
        if (intercept) {
            my <- ifelse(weighted, weighted.mean(y, weights), mean(y))
            coefs <- matrix(my, nrow = 1L, ncol = length(lambda))
            resids <- matrix(y - my, nrow = dX[1L], ncol = length(lambda))
        } else {
            coefs <- matrix(NA_real_, nrow = 0L, ncol = length(lambda))
            resids <- matrix(y, nrow = dX[1L], ncol = length(lambda))
        }
        return(list(
            status = 0L,
            message = "",
            coefficients = coefs,
            residuals = resids
        ))
    }

    elnetres <- if (weighted) {
        .elnet.wfit(
            X,
            y,
            weights,
            alpha,
            lambda,
            intercept,
            addLeading1s,
            options
        )
    } else {
        .elnet.fit(
            X,
            y,
            alpha,
            lambda,
            intercept,
            addLeading1s,
            options
        )
    }

    return(elnetres)
}

## Internal function to fit an EN linear regression WITHOUT parameter checks!
#' @useDynLib pense C_augtrans C_elnet
.elnet.fit <- function(X, y, alpha, lambda, intercept = TRUE,
                       addLeading1s = TRUE, options = en_options_aug_lars(),
                       warm_coefs) {
    y <- drop(y)
    dX <- dim(X)

    ## Add leading column of 1's
    Xtr <- if (!identical(addLeading1s, FALSE)) {
        .Call(C_augtrans, X)
    } else {
        t(X)
    }

    if (missing(warm_coefs)) {
        warm_coefs <- numeric(dX[2L])
        options$warmStart <- FALSE
    } else {
        options$warmStart <- TRUE
    }

    elnetres <- .Call(
        C_elnet,
        Xtr,
        y,
        warm_coefs,
        alpha,
        lambda,
        intercept,
        options
    )

    if (!identical(elnetres[[1L]], 0L)) {
        warning("Elastic Net algorithm had non-zero return status.")
    }

    names(elnetres) <- c("status", "message", "coefficients", "residuals")

    return(elnetres)
}

## Internal function to fit an EN linear regression WITHOUT parameter checks!
#' @useDynLib pense C_augtrans C_elnet_weighted
.elnet.wfit <- function(X, y, weights, alpha, lambda, intercept = TRUE,
                        addLeading1s = TRUE, options = en_options_aug_lars(),
                        warm_coefs) {
    y <- drop(y)
    dX <- dim(X)

    ## Add leading column of 1's
    Xtr <- if (!identical(addLeading1s, FALSE)) {
        .Call(C_augtrans, X)
    } else {
        t(X)
    }

    if (missing(warm_coefs)) {
        warm_coefs <- numeric(dX[2L])
        options$warmStart <- FALSE
    } else {
        options$warmStart <- TRUE
    }

    elnetres <- .Call(
        C_elnet_weighted,
        Xtr,
        y,
        weights,
        warm_coefs,
        alpha,
        lambda,
        intercept,
        options
    )

    if (!identical(elnetres[[1L]], 0L)) {
        warning("Elastic Net algorithm had non-zero return status.")
    }

    names(elnetres) <- c("status", "message", "coefficients", "residuals")

    return(elnetres)
}
