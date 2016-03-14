#' Extract Model Coefficients
#'
#' @param object an object of type \code{pense} to extract coefficients from.
#' @param lambda the value of the penalty parameter. Default is to use the optimal
#'      lambda \code{lambda.opt}.
#' @param ... currently ignored.
#' @return A numeric vector of size \eqn{p + 1}.
#' @export
coef.pense <- function(object, lambda, ...) {
    if (missing(lambda)) {
        return(object$coefficients)
    }

    if (length(lambda) != 1L || !is.numeric(lambda) || is.na(lambda) || lambda < 0) {
        stop("`lambda` must be single non-negative number")
    }

    x <- data.matrix(eval(object$call$x));
    y <- drop(eval(object$call$y));

    ##
    ## Get an initial estimate at the supplied lambda
    ##

    return(NA_real_)
}
