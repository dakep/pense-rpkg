#' Extract Model Coefficients
#'
#' @param object object of type \code{pense} to extract coefficients from.
#' @param lambda the value of the penalty parameter. Default is to use the
#'      optimal lambda \code{lambda.opt}.
#' @param ... currently not used.
#' @return A numeric vector of size \eqn{p + 1}.
#' @export
coef.pense <- function(object, lambda, ...) {
    if (missing(lambda) || is.null(lambda)) {
        return(object$coefficients)
    }

    lambda <- .check_arg(lambda, "numeric", range = 0)

    X <- data.matrix(eval(object$call$X))
    y <- drop(eval(object$call$y))

    estimate <- pense_coldwarm(
        X,
        y,
        lambda_grid = lambda,
        alpha = object$alpha,
        standardize = object$standardize,
        pense_options = object$pense_options,
        initest_options = object$initest_options,
        en_options = object$en_options
    )

    return(nameCoefVec(c(estimate[[1L]]$intercept, estimate[[1L]]$beta), X))
}
