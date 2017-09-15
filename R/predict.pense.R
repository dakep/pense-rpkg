#' Predict Method for Penalized Elastic Net S- and MM-estimators
#'
#' @param object an object of type \code{pense} or \code{pensem} to use for
#'      prediction.
#' @param newdata an optional design matrix
#' @param lambda the value of the penalty parameter. Default is to use the
#'      optimal lambda \code{lambda_opt}.
#' @param exact if the lambda is not part of the lambda grid, should the
#'      estimates be obtained by linear interpolation between the nearest
#'      lambda values (default) or computed exactly.
#' @param correction should a correction factor be applied to the PENSE(M)
#'       estimate?
#' @param ... currently ignored.
#' @return a numeric vector of predicted values for the given lambda.
#' @export
#' @importFrom Matrix drop
predict.pense <- function(object, newdata, lambda, exact = FALSE,
                          correction = TRUE, ...) {
    if (missing(newdata)) {
        newdata <- data.matrix(eval(object$call$X))
    } else {
        dX <- dim(newdata)
        if (is.data.frame(newdata)) {
            newdata <- data.matrix(newdata)
        }

        if (!is.matrix(newdata) || !is.numeric(newdata)) {
            stop("`newdata` must be a numeric matrix")
        }

        if (dX[2L] != nrow(object$coefficients) - 1L) {
            stop("`newdata` must have as many columns as the original data set")
        }
    }

    if (missing(lambda)) {
        lambda <- NULL
    }

    coefs <- coef(object, lambda = lambda, exact = exact, sparse = TRUE,
                  correction = correction)

    pr <- drop(coefs[1L] + newdata %*% coefs[-1L, , drop = FALSE])
    return(as.vector(pr))
}

#' Predict Method for the classical Elastic Net Estimator
#'
#' @param object an object of type \code{elnetfit} to use for prediction.
#' @param newdata an optional design matrix
#' @param lambda the value of the penalty parameter. Default is to use the
#'      optimal lambda \code{lambda_opt}.
#' @param exact if the lambda is not part of the lambda grid, should the
#'      estimates be obtained by linear interpolation between the nearest
#'      lambda values (default) or computed exactly.
#' @param correction should a correction factor be applied to the EN estimate?
#' @param ... currently ignored.
#' @return a numeric vector of predicted values for the given lambda.
#' @export
#' @importFrom Matrix drop
predict.elnetfit <- function(object, newdata, lambda, exact = FALSE,
                             correction = TRUE, ...) {
    if (missing(newdata)) {
        newdata <- data.matrix(eval(object$call$X))
    } else {
        dX <- dim(newdata)
        if (is.data.frame(newdata)) {
            newdata <- data.matrix(newdata)
        }

        if (!is.matrix(newdata) || !is.numeric(newdata)) {
            stop("`newdata` must be a numeric matrix")
        }

        if (dX[2L] != nrow(object$coefficients) - 1L) {
            stop("`newdata` must have as many columns as the original data set")
        }
    }

    if (missing(lambda)) {
        lambda <- NULL
    }

    coefs <- coef.elnetfit(object, lambda = lambda, exact = exact,
                           sparse = TRUE, correction = correction)

    pr <- drop(coefs[1L] + newdata %*% coefs[-1L, , drop = FALSE])
    return(as.vector(pr))
}
