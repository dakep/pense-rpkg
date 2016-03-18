#' Predict Method for Penalized Elastic-Net S-estimator
#'
#' @param object an object of type \code{pense} to use for prediction.
#' @param newdata an optional design matrix
#' @param lambda the value of the penalty parameter. Default is to use the optimal
#'      lambda \code{lambda.opt}.
#' @param ... currently ignored.
#' @return a numeric vector of predicted values for the given lambda.
#' @export
predict.pense <- function(object, newdata, lambda, ...) {
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

        if (dX[2L] != length(object$coefficients) - 1L) {
            stop("`newdata` must have as many columns as the original data set")
        }
    }

    if (missing(lambda)) {
        lambda <- NULL
    }

    coefs <- coef.pense(object, lambda = lambda)

    pr <- coefs[1L] + newdata %*% coefs[-1L]
    return(as.vector(pr))
}
