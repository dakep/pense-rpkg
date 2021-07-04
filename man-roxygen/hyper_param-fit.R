#' @param lambda a single number for the penalty level.
#' @param alpha Either a single number or `NULL` (default).
#'    If given, only fits with the given `alpha` value are considered.
#'    If `object` was fit with multiple `alpha` values, and no value is provided, the
#'    first value in `object$alpha` is used with a warning.
