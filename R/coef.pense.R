#' Extract Model Coefficients
#'
#' @param object an object of type \code{pense} to extract coefficients from.
#' @param lambda the value of the penalty parameter. Default is to use the optimal
#'      lambda \code{lambda.opt}.
#' @param ... currently ignored.
#' @return A numeric vector of size \eqn{p + 1}.
#' @export
coef.pense <- function(object, lambda, ...) {
    if (missing(lambda) || is.null(lambda)) {
        return(object$coefficients)
    }

    if (length(lambda) != 1L || !is.numeric(lambda) || is.na(lambda) || lambda < 0) {
        stop("`lambda` must be single non-negative number")
    }

    X <- data.matrix(eval(object$call$X))
    y <- drop(eval(object$call$y))

    estimate <- pense.coldwarm(X, y, lambda.grid = lambda,
                               alpha = object$alpha,
                               standardize = object$standardize,
                               control = object$control)

    return(nameCoefVec(c(estimate[[1L]]$intercept, estimate[[1L]]$beta), X))
}
