#' Extract Residuals from a Fitted Penalized Elastic-Net S/MM-estimator
#'
#' @param object a PENSE or PENSEM estimate to extract the residuals from.
#' @param lambda the value of the penalty parameter. Default is to use the
#'      optimal lambda \code{object$lambda_opt}.
#' @param exact if the lambda is not part of the lambda grid, should the
#'      estimates be obtained by linear interpolation between the nearest
#'      lambda values (default) or computed exactly.
#' @param correction should a correction factor be applied to the EN estimate?
#'       See \code{\link{elnet}} for details on the applied correction.
#' @param ... currently ignored.
#' @return a numeric vector of residuals for the given lambda.
#' @importFrom Matrix drop
#'
#' @example examples/pense-methods.R
#'
#' @export
residuals.pense <- function(object, lambda, exact = FALSE, correction = TRUE, ...) {
    exact <- isTRUE(exact)

    if (missing(lambda) || is.null(lambda)) {
        lambda <- object$lambda_opt
        exact <- FALSE
    }

    lambda_diff_abs <- abs(object$lambda - lambda)
    lambda_match <- which(lambda_diff_abs < .Machine$double.eps)

    if (length(lambda_match) > 0L) {
        return(object$residuals[ , lambda_match[1L]])
    }

    coefs <- coef.pense(
        object,
        lambda = lambda,
        exact = exact,
        sparse = TRUE,
        correction = correction
    )

    x <- if (!is.null(object$call$x_train)) {
        data.matrix(eval(object$call$x_train))
    } else {
        data.matrix(eval(object$call$x))
    }
    y <- if (!is.null(object$call$y_train)) {
        data.matrix(eval(object$call$y_train))
    } else {
        data.matrix(eval(object$call$y))
    }

    return(drop(y - coefs[1L] - x %*% coefs[-1L, , drop = FALSE]))
}


#' Extract Residuals from a Fitted Elastic-Net Estimator
#'
#' @param object a \code{\link{elnet}} estimate to extract the residuals from.
#' @param lambda the value of the penalty parameter. Default is to use the
#'      optimal lambda \code{object$lambda_opt} if the given estimator was
#'      cross-validated or the smallest lambda if not.
#' @param exact if the lambda is not part of the lambda grid, should the
#'      estimates be obtained by linear interpolation between the nearest
#'      lambda values (default) or computed exactly.
#' @param ... currently ignored.
#' @return a numeric vector of residuals for the given lambda.
#'
#' @importFrom Matrix drop
#'
#' @example examples/elnet_cv-methods.R
#'
#' @export
residuals.elnetfit <- function(object, lambda, exact = FALSE, ...) {
    exact <- isTRUE(exact)

    if (missing(lambda) || is.null(lambda)) {
        if ("cv_elnetfit" %in% class(object)) {
            lambda <- object$lambda_opt
            exact <- FALSE
        } else {
            warning("Using smallest lambda since no CV prediction performance",
                    "is available.")
            lambda <- object$lambda[1L]
            exact <- FALSE
        }
    }

    lambda_diff_abs <- abs(object$lambda - lambda)
    lambda_match <- which(lambda_diff_abs < .Machine$double.eps)

    if (length(lambda_match) > 0L) {
        return(object$residuals[ , lambda_match[1L]])
    }

    coefs <- coef.elnetfit(
        object,
        lambda = lambda,
        exact = exact,
        sparse = TRUE
    )

    x <- data.matrix(eval(object$call$x))
    y <- drop(eval(object$call$y))

    return(drop(y - coefs[1L] - x %*% coefs[-1L, , drop = FALSE]))
}


