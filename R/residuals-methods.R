#' Predict Method for PENSE Fits
#'
#' Predict response values using a PENSE (or LS-EN) regularization path fitted by [pense()] or [elnet()].
#'
#' @param object PENSE regularization path to extract residuals from.
#' @param newdata an optional matrix of new predictor values. If missing, the fitted values are computed.
#' @param lambda a single value of the penalty parameter.
#' @param exact defunct Always gives a warning if `lambda` is not part of the fitted sequence and coefficients
#'    need to be interpolated.
#' @param correction defunct.
#' @param ... currently not used.
#' @return a numeric vector of residuals for the given penalization level.
#'
#' @family functions for extracting components
#'
#' @example examples/residuals.R
#' @export
#'
#' @importFrom Matrix drop
#' @importFrom lifecycle deprecate_warn deprecated is_present
#' @importFrom rlang abort
#' @importFrom stats coef
#'
predict.pense_fit <- function(object, newdata, lambda, exact = deprecated(), correction = deprecated(), ...) {
  if (is_present(exact)) {
    deprecate_warn('2.0.0', 'residuals(exact=)')
  }
  if (is_present(correction)) {
    deprecate_stop('2.0.0', 'residuals(correction=)')
  }

  coef_est <- coef(object, lambda = lambda, concat = FALSE)

  if (missing(newdata)) {
    newdata <- eval.parent(object$call$x)
  } else {
    newdata <- data.matrix(newdata)
    x_dim <- dim(newdata)
    if (x_dim[[2L]] != length(coef_est$beta)) {
      abort(sprintf("`newdata` is of invalid size (contains %d predictors instead of the required %d)",
                    x_dim[[2L]], length(coef_est$beta)))
    }
  }

  return(drop(newdata %*% coef_est$beta) + coef_est$intercept)
}

#' Predict Method for PENSE Fits
#'
#' Predict response values using a PENSE (or LS-EN) regularization path with hyper-parameters chosen by
#' cross-validation.
#'
#' If `lambda = "se"` and `object` contains fitted estimates for every penalization level in the sequence, extract the
#' residuals of the most parsimonious model with prediction performance statistically indistinguishable from the best
#' model. This is determined to be the model with prediction performance within `se_mult * cv_se` from the best model.
#'
#' @param object PENSE with cross-validated hyper-parameters to extract coefficients from.
#' @param newdata an optional matrix of new predictor values. If missing, the fitted values are computed.
#' @param lambda either a string specifying which penalty level to use or a a single numeric value of the penalty
#'    parameter. See details.
#' @param se_mult If `lambda = "se"`, the multiple of standard errors to tolerate.
#' @param exact deprecated. Always gives a warning if `lambda` is not part of the fitted sequence and coefficients
#'    are interpolated.
#' @param correction defunct.
#' @param ... currently not used.
#' @return a numeric vector of residuals for the given penalization level.
#'
#' @family functions for extracting components
#'
#' @example examples/residuals.R
#' @export
#'
#' @importFrom Matrix drop
#' @importFrom lifecycle deprecate_warn deprecated is_present
#' @importFrom rlang warn
#' @importFrom stats coef
#'
predict.pense_cvfit <- function(object, newdata, lambda = c('min', 'se'), se_mult = 1, exact = deprecated(),
                                correction = deprecated(), ...) {
  if (is_present(exact)) {
    deprecate_warn('2.0.0', 'coef(exact=)')
  }
  if (is_present(correction)) {
    deprecate_stop('2.0.0', 'coef(correction=)')
  }
  if (is.character(lambda)) {
    lambda <- match.arg(lambda)
  }

  coef_est <- coef(object, lambda = lambda, se_mult = se_mult, concat = FALSE)

  if (missing(newdata)) {
    newdata <- eval.parent(object$call$x)
  } else {
    newdata <- data.matrix(newdata)
    x_dim <- dim(newdata)
    if (x_dim[[2L]] != length(coef_est$beta)) {
      abort(sprintf("`newdata` is of invalid size (contains %d predictors instead of the required %d)",
                    x_dim[[2L]], length(coef_est$beta)))
    }
  }

  return(drop(newdata %*% coef_est$beta) + coef_est$intercept)
}

#' Extract Residuals
#'
#' Extract residuals from a PENSE (or LS-EN) regularization path fitted by [pense()] or [elnet()].
#'
#' @param object PENSE regularization path to extract residuals from.
#' @param lambda a single value of the penalty parameter.
#' @param exact defunct Always gives a warning if `lambda` is not part of the fitted sequence and coefficients
#'    need to be interpolated.
#' @param correction defunct.
#' @param ... currently not used.
#' @return a numeric vector of residuals for the given penalization level.
#'
#' @family functions for extracting components
#'
#' @example examples/residuals.R
#' @export
residuals.pense_fit <- function(object, lambda, exact = deprecated(), correction = deprecated(), ...) {
  train_y <- eval.parent(object$call$y)
  cl <- match.call(expand.dots = FALSE)
  cl[[1L]] <- quote(predict)
  return(train_y - eval.parent(cl))
}

#' Extract Residuals
#'
#' Extract residuals from a PENSE (or LS-EN) regularization path with hyper-parameters chosen by cross-validation.
#'
#' If `lambda = "se"` and `object` contains fitted estimates for every penalization level in the sequence, extract the
#' residuals of the most parsimonious model with prediction performance statistically indistinguishable from the best
#' model. This is determined to be the model with prediction performance within `se_mult * cv_se` from the best model.
#'
#' @param object PENSE with cross-validated hyper-parameters to extract coefficients from.
#' @param lambda either a string specifying which penalty level to use or a a single numeric value of the penalty
#'    parameter. See details.
#' @param se_mult If `lambda = "se"`, the multiple of standard errors to tolerate.
#' @param exact deprecated. Always gives a warning if `lambda` is not part of the fitted sequence and coefficients
#'    are interpolated.
#' @param correction defunct.
#' @param ... currently not used.
#' @return a numeric vector of residuals for the given penalization level.
#'
#' @family functions for extracting components
#'
#' @example examples/residuals.R
#' @export
residuals.pense_cvfit <- function(object, lambda = c('min', 'se'), se_mult = 1, exact = deprecated(),
                                  correction = deprecated(), ...) {
  train_y <- eval.parent(object$call$y)
  cl <- match.call(expand.dots = FALSE)
  cl[[1L]] <- quote(predict)
  return(train_y - eval.parent(cl))
}
