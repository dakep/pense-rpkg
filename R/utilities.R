#' Compute the Tau-Scale of Centered Values
#'
#' Compute the \eqn{\tau}-scale without centering the values.
#'
#' @param x numeric values. Missing values are verbosely ignored.
#' @return the \eqn{\tau} estimate of scale of centered values.
#'
#' @family functions to compute robust estimates of location and scale
#'
#' @export
#' @importFrom rlang warn
#' @importFrom stats na.omit
tau_size <- function (x) {
  x <- if (anyNA(x)) {
    warn("Missing values are ignored.")
    .as(na.omit(x), 'numeric')
  } else {
    .as(x, 'numeric')
  }
  .Call(C_tau_size, x)
}

#' Compute the M-Scale of Centered Values
#'
#' Compute the M-scale without centering the values.
#'
#' @param x numeric values. Missing values are verbosely ignored.
#' @param bdp desired breakdown point (between 0 and 0.5).
#' @param cc cutoff value for the bisquare rho function.
#'    By default, chosen to yield a consistent estimate for the Normal distribution.
#' @param opts a list of options for the M-scale estimation algorithm,
#'    see [mscale_algorithm_options()] for details.
#' @param delta deprecated. Use `bpd` instead.
#' @param rho,eps,maxit deprecated. Instead set control options for the algorithm
#'    with the `opts` arguments.
#' @return the M-estimate of scale.
#'
#' @family functions to compute robust estimates of location and scale
#'
#' @export
#'
#' @importFrom lifecycle deprecate_stop deprecated is_present
#' @importFrom rlang warn
#' @importFrom stats na.omit
mscale <- function (x, bdp = 0.25, cc = consistency_const(bdp, 'bisquare'),
                    opts = mscale_algorithm_options(), delta = deprecated(),
                    rho = deprecated(), eps = deprecated(), maxit = deprecated()) {
  if (is_present(delta)) {
    deprecate_stop('2.0.0', 'mscale(delta=)', 'mscale(bdp=)')
    bdp <- delta
  }
  if (is_present(rho)) {
    deprecate_stop('2.0.0', 'mscale(rho=)', 'mscale(opts=)')
  }
  if (is_present(eps)) {
    deprecate_stop('2.0.0', 'mscale(rho=)', 'mscale(opts=)')
    opts$eps <- .as(eps[[1L]], 'numeric')
  }
  if (is_present(maxit)) {
    deprecate_stop('2.0.0', 'mscale(maxit=)', 'mscale(opts=)')
    opts$max_it <- .as(maxit[[1L]], 'integer')
  }

  x <- if (anyNA(x)) {
    warn("Missing values are ignored.")
    .as(na.omit(x), 'numeric')
  } else {
    .as(x, 'numeric')
  }

  if (missing(cc)) {
    cc <- NULL
  }
  opts <- .full_mscale_algo_options(bdp, cc, opts)
  .Call(C_mscale, x, opts)
}

#' Compute the Gradient and Hessian of the M-Scale Function
#'
#' Compute the derivative (gradient) or the Hessian of the M-scale function
#' evaluated at the point `x`.
#'
#' @param x numeric values. Missing values are verbosely ignored.
#' @param bdp desired breakdown point (between 0 and 0.5).
#' @param cc cutoff value for the bisquare rho function.
#'    By default, chosen to yield a consistent estimate for the
#'    Normal distribution.
#' @param order compute the gradient (`order=1`) or the gradient and the
#'    Hessian (`order=2`).
#' @param opts a list of options for the M-scale estimation algorithm,
#'    see [mscale_algorithm_options()] for details.
#' @return a vector of derivatives of the M-scale function, one per element in `x`.
#'
#' @importFrom rlang warn
#' @importFrom stats na.omit
#' @keywords internal
mscale_derivative <- function (x, bdp = 0.25, order = 1,
                               cc = consistency_const(bdp, 'bisquare'),
                               opts = mscale_algorithm_options()) {
  x <- if (anyNA(x)) {
    warn("Missing values are ignored.")
    .as(na.omit(x), 'numeric')
  } else {
    .as(x, 'numeric')
  }

  order <- .as(order[[1L]], 'integer')
  if (order < 1L) {
    order <- 1L
  } else if (order > 2L) {
    order <- 2L
  }

  if (missing(cc)) {
    cc <- NULL
  }
  opts <- .full_mscale_algo_options(bdp, cc, opts)
  grad_hess <- .Call(C_mscale_derivative, x, opts, order)

  if (identical(order, 1L)) {
    if (length(grad_hess) < length(x)) {
      grad_hess <- rep_len(NA_real_, length(x))
    }
    return(grad_hess)
  } else {
    if (length(grad_hess) < length(x)) {
      return(list(gradient = rep_len(NA_real_, length(x)),
                  rho_2nd = rep_len(NA_real_, length(x)),
                  denom = NA_real_,
                  scale = grad_hess[[1L]],
                  violation = NA_real_,
                  hessian = matrix(NA_real_, length(x), length(x))))
    }
    grad_hess <- list(gradient = grad_hess[, 1L],
                      rho_2nd = grad_hess[, 2L],
                      denom = grad_hess[2L, 3L],
                      scale = grad_hess[3L, 3L],
                      violation = grad_hess[4L, 3L],
                      hessian = grad_hess[, -(1:2)])
    grad_hess$hessian[2:4, 1L] <- 0
    return(grad_hess)
  }
}

#' @description
#' Compute the maximum derivative of the M-scale function with respect to each element over
#' a grid of values.
#'
#' @param n_change the number of elements in `x` to replace with each value in `grid`.
#' @param grid a grid of values to replace the first 1 - `n_change` elements in` x`.
#' @return a vector with 4 elements:
#'    1. the maximum absolute value of the gradient,
#'    2. the maximum absolute value of the Hessian elements,
#'    3. the M-scale associated with 1., and
#'    4. the M-scale associated with 2.
#' @describeIn mscale_derivative maximum of the gradient
#' @keywords internal
max_mscale_derivative <- function (x, grid, n_change, bdp = 0.25,
                                   cc = consistency_const(bdp, 'bisquare'),
                                   opts = mscale_algorithm_options()) {
  x <- if (anyNA(x)) {
    warn("Missing values are ignored.")
    .as(na.omit(x), 'numeric')
  } else {
    .as(x, 'numeric')
  }

  if (missing(cc)) {
    cc <- NULL
  }
  opts <- .full_mscale_algo_options(bdp, cc, opts)

  .Call(C_max_mscale_derivative, x, grid, n_change, opts)
}

#' @description
#' Compute the maximum element in the gradient and Hessian of the M-scale
#' function with respect to each element over a grid of values.
#'
#' @param n_change the number of elements in `x` to replace with each value in `grid`.
#' @param grid a grid of values to replace the first 1 - `n_change` elements in` x`.
#' @return the maximum absolute derivative over the entire grid.
#' @describeIn mscale_derivative maximum of the gradient and hessian
#' @keywords internal
max_mscale_grad_hess <- function (x, grid, n_change, bdp = 0.25,
                                  cc = consistency_const(bdp, 'bisquare'),
                                  opts = mscale_algorithm_options()) {
  x <- if (anyNA(x)) {
    warn("Missing values are ignored.")
    .as(na.omit(x), 'numeric')
  } else {
    .as(x, 'numeric')
  }

  if (missing(cc)) {
    cc <- NULL
  }
  opts <- .full_mscale_algo_options(bdp, cc, opts)

  .Call(C_max_mscale_grad_hess, x, grid, n_change, opts)
}

#' Compute the M-estimate of Location
#'
#' Compute the M-estimate of location using an auxiliary estimate of the scale.
#'
#' @param x numeric values. Missing values are verbosely ignored.
#' @param scale scale of the `x` values. If omitted, uses the [mad()][stats::mad()].
#' @param rho the \eqn{\rho} function to use. See [rho_function()] for available functions.
#' @param cc value of the tuning constant for the chosen \eqn{\rho} function.
#'    By default, chosen to achieve 95% efficiency under the Normal distribution.
#' @param opts a list of options for the M-estimating algorithm, see
#'    [mscale_algorithm_options()] for details.
#' @return a single numeric value, the M-estimate of location.
#'
#' @family functions to compute robust estimates of location and scale
#'
#' @export
#'
#' @importFrom stats mad
#' @importFrom rlang warn
#' @importFrom stats na.omit
mloc <- function (x, scale, rho, cc, opts = mscale_algorithm_options()) {
  x <- if (anyNA(x)) {
    warn("Missing values are ignored.")
    .as(na.omit(x), 'numeric')
  } else {
    .as(x, 'numeric')
  }
  if (missing(scale)) {
    scale <- mad(x)
  }

  if (!isTRUE(scale > .Machine$double.eps)) {
    warn("Cannot compute M-estimate of location for values with scale of 0.")
    return(NA_real_)
  }

  if (missing(cc)) {
    cc <- NULL
  }
  opts <- .full_mscale_algo_options(.5, cc, opts)
  opts$rho <- rho_function(rho)
  .Call(C_mloc, .as(x, 'numeric'), scale, opts)
}

#' Compute the M-estimate of Location and Scale
#'
#' Simultaneous estimation of the location and scale by means of M-estimates.
#'
#' @param x numeric values. Missing values are verbosely ignored.
#' @param bdp desired breakdown point (between 0 and 0.5).
#' @param scale_cc cutoff value for the bisquare \eqn{\rho} function for computing the
#'    scale estimate.
#'    By default, chosen to yield a consistent estimate for normally distributed values.
#' @param location_rho,location_cc \eqn{\rho} function and cutoff value for computing
#'    the location estimate.
#'    See [rho_function()] for a list of available \eqn{\rho} functions.
#' @param opts a list of options for the M-estimating equation,
#'    see [mscale_algorithm_options()] for details.
#' @return a vector with 2 elements, the M-estimate of location and the M-scale estimate.
#'
#' @family functions to compute robust estimates of location and scale
#'
#' @export
#'
#' @importFrom rlang warn
#' @importFrom stats na.omit
mlocscale <- function (x, bdp = 0.25, scale_cc = consistency_const(bdp, 'bisquare'), location_rho,
                       location_cc, opts = mscale_algorithm_options()) {
  x <- if (anyNA(x)) {
    warn("Missing values are ignored.")
    .as(na.omit(x), 'numeric')
  } else {
    .as(x, 'numeric')
  }

  opts <- .full_mscale_algo_options(bdp, scale_cc, opts)
  loc_opts <- list(rho = rho_function(location_rho))
  if (!missing(location_cc)) {
    loc_opts$cc <- .as(location_cc[[1L]], 'numeric')
  }
  .Call(C_mlocscale, x, opts, loc_opts)
}

#' Get the Constant for Consistency for the M-Scale
#'
#' @param delta desired breakdown point (between 0 and 0.5)
#' @param rho the name of the chosen \eqn{\rho} function.
#'
#' @return consistency constant
#'
#' @family miscellaneous functions
#'
#' @export
#'
#' @importFrom rlang abort
consistency_const <- function (delta, rho) {
  return(switch(rho_function(rho),
                bisquare = .bisquare_consistency_const(delta),
                huber = abort("Huber's rho function not supported for scale estimation!")))
}

#' List Available Rho Functions
#'
#' @param rho the name of the \eqn{\rho} function to check for existence.
#' @return if `rho` is missing returns a vector of supported \eqn{\rho} function names, otherwise
#'    the internal integer representation of the \eqn{\rho} function.
#'
#' @family miscellaneous functions
#'
#' @export
rho_function <- function (rho) {
  available <- c('bisquare', 'huber')
  if (missing(rho)) {
    return(available)
  }
  return(match(match.arg(rho, available), available))
}

#' Create Starting Points for the PENSE Algorithm
#'
#' Create a starting point for starting the PENSE algorithm in [pense()].
#' Multiple starting points can be created by combining starting points via
#' `c(starting_point_1, starting_point_2, ...)`.
#'
#' A starting points can either be *shared*, i.e., used for every penalization level PENSE
#' estimates are computed for, or *specific* to one penalization level.
#' To create a specific starting point, provide the penalization parameters `lambda` and `alpha`.
#' If `lambda` or `alpha` are missing, a shared starting point is created.
#' Shared and specific starting points can all be combined into a single list of starting points,
#' with [pense()] handling them correctly.
#' Note that specific starting points will lead to the `lambda` value being added to the
#' grid of penalization levels.
#' See [pense()] for details.
#'
#' Starting points computed via [enpy_initial_estimates()] are by default *shared* starting points
#' but can be transformed to *specific* starting points via
#' `as_starting_point(..., specific = TRUE)`.
#'
#' @param beta beta coefficients at the starting point. Can be a numeric vector,
#'    a sparse vector of class [dsparseVector][Matrix::sparseVector-class],
#'    or a sparse matrix of class [dgCMatrix][Matrix::CsparseMatrix-class] with a single column.
#' @param intercept intercept coefficient at the starting point.
#' @param alpha optional value(s) for the `alpha` hyper-parameter.
#'    For `starting_point()`, makes the starting point specific to the given `alpha` value.
#'    For `as_starting_point()`, only estimates with matching `alpha` values are extracted.
#'    Values not present in `object` are ignored with a warning.
#' @param lambda optional penalization level(s).
#'    For `starting_point()`, makes the starting point specific to the given `alpha` value.
#'    For `as_starting_point()`, only estimates with matching `alpha` and `lambda` values
#'    are extracted.
#'    If `lambda` is given, `alpha` must be given and be a single number.
#'    Penalization levels not present in `object` are ignored with a warning.
#' @return an object of type `starting_points` to be used as starting point for [pense()].
#'
#' @family functions for initial estimates
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom Matrix sparseVector
#' @importFrom rlang abort
starting_point <- function (beta, intercept, lambda, alpha) {
  sp <- list(intercept = .as(intercept[[1L]], 'numeric'))
  if (missing(lambda) && missing(alpha)) {
    class(sp) <- c('shared_starting_point', 'starting_point')
  } else {
    sp$lambda <- .as(lambda[[1L]], 'numeric')
    sp$alpha <- .as(alpha[[1L]], 'numeric')
    class(sp) <- c('specific_starting_point', 'starting_point')
  }
  if (is(beta, 'dgCMatrix') && ncol(beta) == 1L) {
    beta <- sparseVector(beta@x, beta@i + 1L, beta@Dim[[1L]])
  }

  if (is(beta, 'dsparseVector') || is(beta, 'numeric')) {
    sp$beta <- beta
  } else {
    abort("`beta` must be a (sparse) numeric vector.")
  }
  return(sp)
}

#' @rdname starting_point
#' @param object an object with estimates to use as starting points.
#' @param specific whether the estimates should be used as starting points only at
#'   the penalization level they are computed for.
#'   Defaults to using the estimates as starting points for all penalization levels.
#' @param ... further arguments passed to or from other methods.
#' @export
as_starting_point <- function (object, specific = FALSE, ...) {
  UseMethod('as_starting_point')
}

#' @rdname starting_point
#' @export
as_starting_point.enpy_starting_points <- function (object, specific = FALSE, ...) {
  if (isTRUE(specific)) {
    structure(lapply(object, structure, class = c('specific_starting_point', 'starting_point')),
              class = 'starting_points')
  } else {
    object
  }
}

#' @rdname starting_point
#' @export
#' @importFrom rlang warn abort
as_starting_point.pense_fit <- function (object, specific = FALSE, alpha, lambda, ...) {
  alpha <- if (missing(alpha) || is.null(alpha)) {
    orig_alpha_length <- length(object$alpha)
    object$alpha
  } else {
    orig_alpha_length <- length(alpha)
    object$alpha[na.omit(.approx_match(.as(alpha, 'numeric'), object$alpha))]
  }
  if (length(alpha) == 0L) {
    abort("No requested `alpha` values are not available in `object`.")
  } else if (length(alpha) < orig_alpha_length) {
    warn("Some requested `alpha` values are not available in `object`.")
  }

  alpha_indices <- which(vapply(object$estimates, FUN.VALUE = logical(1L), FUN = function (est) {
    any((est$alpha - alpha)^2 < .Machine$double.eps)
  }))

  lambda_indices <- if (missing(lambda)) {
    alpha_indices
  } else if (is.numeric(lambda)) {
    if (length(alpha) > 1L) {
      abort("If `lambda` is given `alpha` must be a single number.")
    }
    ai_lambdas <- vapply(object$estimates[alpha_indices], FUN.VALUE = numeric(1L),
                         FUN = `[[`, 'lambda')
    lambda_indices <- .approx_match(lambda, ai_lambdas)
    bad_lambda_indices <- which(is.na(lambda_indices))
    if (length(bad_lambda_indices) == length(lambda)) {
      abort("No penalization level requested with `lambda` is available in `object`.")
    } else if (length(bad_lambda_indices) > 0L) {
      warn("Some penalization levels requested with `lambda` are not available in `object`.")
      alpha_indices[lambda_indices[-bad_lambda_indices]]
    } else {
      alpha_indices[lambda_indices]
    }
  } else {
    abort("`lambda` must be missing or a numeric vector.")
  }

  class_list <- if (isTRUE(specific)) {
    c('specific_starting_point', 'starting_point')
  } else {
    c('shared_starting_point', 'starting_point')
  }

  structure(lapply(object$estimates[lambda_indices], structure, class = class_list),
            class = 'starting_points')
}

#' @rdname starting_point
#' @param alpha optional value for the `alpha` hyper-parameter.
#'    If given, only estimates with matching `alpha` values are extracted.
#'    Values not present in `object` are ignored with a warning.
#' @param lambda optionally either a string specifying which penalty level to use
#'    (`"min"` or `"se"`) or a numeric vector of the penalty levels to extract from `object`.
#'    Penalization levels not present in `object` are ignored with a warning.
#'    If `NULL`, all estimates in `object` are extracted.
#'    If a numeric vector, `alpha` must be given and a single number.
#' @param se_mult If `lambda = "se"`, the multiple of standard errors to tolerate.
#'
#' @details
#' When creating starting points from cross-validated fits, it is possible to extract only the
#' estimate with best CV performance (`lambda = "min"`), or the estimate with CV performance
#' statistically indistinguishable from the best performance (`lambda = "se"`).
#' This is determined to be the estimate with prediction performance within
#' `se_mult * cv_se` from the best model.
#'
#' @export
#' @importFrom rlang warn abort
as_starting_point.pense_cvfit <- function (object, specific = FALSE,
                                           alpha, lambda = c('min', 'se'),
                                           se_mult = 1, ...) {
  alpha <- if (missing(alpha) || is.null(alpha)) {
    orig_alpha_length <- length(object$alpha)
    object$alpha
  } else {
    orig_alpha_length <- length(alpha)
    object$alpha[na.omit(.approx_match(.as(alpha, 'numeric'), object$alpha))]
  }
  if (length(alpha) == 0L) {
    abort("No requested `alpha` values are not available in `object`.")
  } else if (length(alpha) < orig_alpha_length) {
    warn("Some requested `alpha` values are not available in `object`.")
  }

  alpha_indices <- which(vapply(object$estimates, FUN.VALUE = logical(1L), FUN = function (est) {
    any((est$alpha - alpha)^2 < .Machine$double.eps)
  }))

  lambda_indices <- if (isFALSE(object$call$fit_all)) {
    if (is.character(lambda)) {
      lambda <- match.arg(lambda)
    }
    if (!missing(lambda) && !isTRUE(lambda == 'min')) {
      warn(paste("`object` was created with `fit_all = FALSE`.",
                 "Only the estimates at the minimum is available and",
                 "will be used as starting points."))
    }
    alpha_indices
  } else if (is.character(lambda)) {
    lambda <- match.arg(lambda)
    se_mult <- if (lambda == 'min') {
      0
    } else {
      .as(se_mult[[1L]], 'numeric')
    }

    vapply(alpha, FUN.VALUE = integer(1L), FUN = function (al) {
      rows <- which((object$cvres$alpha - al)^2 < .Machine$double.eps)
      se_selection <- .cv_se_selection(object$cvres$cvavg[rows], object$cvres$cvse[rows], se_mult)
      rows[which(se_selection == 'se_fact')]
    })
  } else if (is.numeric(lambda)) {
    if (length(alpha) > 1L) {
      abort("If `lambda` is numeric, `alpha` must be a single number.")
    }
    ai_lambdas <- vapply(object$estimates[alpha_indices], FUN.VALUE = numeric(1L),
                         FUN = `[[`, 'lambda')
    lambda_indices <- .approx_match(lambda, ai_lambdas)
    bad_lambda_indices <- which(is.na(lambda_indices))
    if (length(bad_lambda_indices) == length(lambda)) {
      abort("No penalization level requested with `lambda` is available in `object`.")
    } else if (length(bad_lambda_indices) > 0L) {
      warn("Some penalization levels requested with `lambda` are not available in `object`.")
      alpha_indices[lambda_indices[-bad_lambda_indices]]
    } else {
      alpha_indices[lambda_indices]
    }
  } else if (is.null(lambda)) {
    alpha_indices
  } else {
    abort("`lambda` must be missing, a numeric vector, or either one of \"min\" or \"se\".")
  }

  class_list <- if (isTRUE(specific)) {
    c('specific_starting_point', 'starting_point')
  } else {
    c('shared_starting_point', 'starting_point')
  }

  structure(lapply(object$estimates[lambda_indices], structure, class = class_list),
            class = 'starting_points')
}

#' @export
#' @keywords internal
#' @importFrom rlang abort
as_starting_point.starting_point <- function (object, specific = FALSE, ...) {
  if (isTRUE(specific)) {
    if (!is.null(object$lambda) && !is.null(object$alpha)) {
      class(object) <- c('specific_starting_point', 'starting_point')
    } else {
      abort(paste("`object` can not be made a specific starting point because it does",
                  "not contain items `lambda` and `alpha`."))
    }
  } else {
    class(object) <- c('shared_starting_point', 'starting_point')
  }
  return(object)
}

#' @export
#' @keywords internal
c.starting_point <- function (...) {
  c.starting_points(...)
}

#' @export
#' @keywords internal
#' @importFrom methods is
#' @importFrom rlang abort
c.starting_points <- function (...) {
  structure(unlist(lapply(list(...), function (sp) {
    if (is(sp, 'starting_point')) {
      list(sp)
    } else if (is(sp, 'starting_points')) {
      sp
    } else {
      abort(sprintf("Do not know how to combine starting points with an object of type `%s`",
                    class(sp)))
    }
  }), recursive = FALSE, use.names = FALSE), class = 'starting_points')
}
