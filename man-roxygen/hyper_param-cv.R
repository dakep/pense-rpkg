#' @section Hyper-parameters:
#'
#' If `lambda = "{m}-se"` and `object` contains fitted estimates for every penalization
#' level in the sequence, use the fit the most parsimonious model with prediction performance
#' statistically indistinguishable from the best model.
#' This is determined to be the model with prediction performance within `m * cv_se`
#' from the best model.
#' If `lambda = "se"`, the multiplier *m* is taken from `se_mult`.
#'
#' By default all *alpha* hyper-parameters available in the fitted object are considered.
#' This can be overridden by supplying one or multiple values in parameter `alpha`.
#' For example, if `lambda = "1-se"` and `alpha` contains two values, the "1-SE" rule is applied
#' individually for each `alpha` value, and the fit with the better prediction error is considered.
#'
#' In case `lambda` is a number and `object` was fit for several *alpha* hyper-parameters,
#' `alpha` must also be given, or the first value in `object$alpha` is used with a warning.
#'
#'
#' @param lambda either a string specifying which penalty level to use
#'    (`"min"`, `"se"`, `"{m}-se`")
#'    or a single numeric value of the penalty parameter. See details.
#' @param alpha Either a single number or `NULL` (default).
#'    If given, only fits with the given `alpha` value are considered.
#'    If `lambda` is a numeric value and `object` was fit with multiple *alpha*
#'    values and no value is provided, the first value in `object$alpha` is used with a warning.
#' @param se_mult If `lambda = "se"`, the multiple of standard errors to tolerate.
