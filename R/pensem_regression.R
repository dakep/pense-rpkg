#' Compute Penalized Elastic Net M-Estimates from PENSE
#'
#' This is a convenience wrapper around [pense_cv()] and [regmest_cv()], for the common use-case of computing
#' a highly-robust S-estimate followed by a more efficient M-estimate using the scale of the residuals from the
#' S-estimate.
#'
#' @param x either a numeric matrix of predictor values, or a cross-validated PENSE fit from [pense_cv()].
#' @family functions to compute robust estimates with CV
#' @export
pensem_cv <- function (x, ...) {
  UseMethod('pensem_cv')
}

#' Deprecated Alias of pensem_cv
#'
#' `pensem()` is a deprecated alias for [pensem_cv()].
#'
#' @inheritParams pensem_cv
#'
#' @family deprecated functions
#'
#' @export
pensem <- function (x, ...) {
  deprecate_warn('2.0.0', 'pensem()', 'pensem_cv()')
  UseMethod('pensem_cv')
}


#' @inheritParams pense
#' @param cc cutoff constant for Tukey's bisquare \eqn{\rho} function in the M-estimation objective function.
#' @param lambda_m,lambda_s optional user-supplied sequence of penalization levels for the S- and M-estimates.
#'    If given and not `NULL`, `nlambda` and `lambda_min_ratio` are ignored for the respective estimate (S and/or M).
#' @template cv_params
#'
#' @return an object of cross-validated regularized M-estimates as returned from [regmest_cv()].
#'
#' @export
#'
#' @rdname pensem_cv
#' @importFrom rlang abort
#' @importFrom methods is
#' @importFrom stats coef
pensem_cv.default <- function (x, y, alpha = 0.5, nlambda = 50, lambda_min_ratio, lambda_m, lambda_s,
                               standardize = TRUE, penalty_loadings, intercept = TRUE, bdp = 0.25, ncores = 1,
                               sparse = FALSE, eps = 1e-6, cc = 4.7, cv_k = 5, cv_repl = 1, cl = NULL,
                               cv_metric = c('tau_size', 'mape', 'rmspe'), add_zero_based = TRUE,
                               explore_solutions = 10, explore_tol = 0.1, max_solutions = 10,
                               fit_all = TRUE, comparison_tol = sqrt(eps), algorithm_opts = mm_algorithm_options(),
                               mscale_opts = mscale_algorithm_options(), nlambda_enpy = 10, enpy_opts = enpy_options(),
                               ...) {
  if (is(x, 'pense_fit')) {
    abort("A `pense()` fit can not be used for `pensem()`. Use `pense_cv()` instead.")
  }

  cl <- match.call(expand.dots = TRUE)
  if (missing(cv_k)) {
    cl$cv_k <- cv_k
  }

  ## Compute the PENSE fit
  pense_call <- cl
  pense_call[[1L]] <- quote(pense::pense_cv)
  pense_call$fit_all <- FALSE
  pense_call$lambda_m <- NULL
  if (!missing(lambda_s)) {
    pense_call$lambda <- lambda_s
    pense_call$lambda_s <- NULL
  }
  pense_fit <- eval.parent(pense_call)

  ## Compute the PENSEM fit
  # Determine the scale of the residuals from the best PENSE fit
  pense_coefs <- coef(pense_fit, lambda = 'min')
  resids <- as.numeric(y - x %*% pense_coefs[-1L] - pense_coefs[[1L]])
  resid_scale <- mscale(resids, bdp = bdp, opts = mscale_opts)

  regmest_call <- cl
  regmest_call[[1L]] <- quote(pense::regmest_cv)
  regmest_call$fit_all <- fit_all
  regmest_call$starting_points <- as_starting_point(pense_fit, lambda = 'min')
  regmest_call$scale <- resid_scale
  regmest_call$mscale_bdp <- bdp
  regmest_call$lambda_s <- NULL
  if (!missing(lambda_m)) {
    regmest_call$lambda <- lambda_m
    regmest_call$lambda_m <- NULL
  }
  regmest_fit <- eval.parent(regmest_call)

  regmest_fit$initial <- pense_fit

  return(regmest_fit)
}

#' @rdname pensem_cv
#'
#' @param scale initial scale estimate to use in the M-estimation. By default the S-scale from the PENSE fit is used.
#' @param x_train,y_train override arguments `x` and `y` as provided in the call to `pense_cv()`. This is useful if
#'    the arguments in the `pense_cv()` call are not available in the current environment.
#'
#' @importFrom stats coef
#' @export
#'
#' @seealso [pense_cv()] to compute the starting S-estimate.
pensem_cv.pense_cvfit <- function (x, scale, alpha, nlambda = 50, lambda_min_ratio, lambda_m,
                                   standardize = TRUE, penalty_loadings, intercept = TRUE, bdp = 0.25, ncores = 1,
                                   sparse = FALSE, eps = 1e-6, cc = 4.7, cv_k = 5, cv_repl = 1, cl = NULL,
                                   cv_metric = c('tau_size', 'mape', 'rmspe'), add_zero_based = TRUE,
                                   explore_solutions = 10, explore_tol = 0.1, max_solutions = 10,
                                   fit_all = TRUE, comparison_tol = sqrt(eps), algorithm_opts = mm_algorithm_options(),
                                   mscale_opts = mscale_algorithm_options(), x_train, y_train, ...) {
  cl <- match.call(expand.dots = TRUE)

  if (missing(cv_k)) {
    cl$cv_k <- cv_k
  }

  # Extract information from pense fit
  if (missing(alpha)) {
    alpha <- x$estimates[[1L]]$alpha
  }

  if (missing(x_train)) {
    x_train <- eval.parent(x$call$x)
  }
  if (missing(y_train)) {
    y_train <- eval.parent(x$call$y)
  }

  ## Compute the PENSEM fit
  # Determine the scale of the residuals from the best PENSE fit
  if (missing(scale)) {
    pense_coefs <- coef(x, lambda = 'min')
    resids <- as.numeric(y_train - x_train %*% pense_coefs[-1L] - pense_coefs[[1L]])
    cc <- if (!is.null(x$call$cc)) {
      eval.parent(x$call$cc)
    } else {
      NULL
    }
    scale <- mscale(resids, bdp = bdp, cc = cc, opts = mscale_opts)
  }

  regmest_call <- cl
  regmest_call[[1L]] <- quote(pense::regmest_cv)
  regmest_call$x <- x_train
  regmest_call$y <- y_train
  regmest_call$alpha <- alpha
  regmest_call$starting_points <- as_starting_point(x, lambda = 'min')
  regmest_call$scale <- scale
  regmest_call$mscale_bdp <- bdp
  regmest_call$lambda_s <- NULL
  regmest_call$x_train <- NULL
  regmest_call$y_train <- NULL
  if (!missing(lambda_m)) {
    regmest_call$lambda <- lambda_m
    regmest_call$lambda_m <- NULL
  }
  regmest_fit <- eval.parent(regmest_call)

  return(regmest_fit)
}
