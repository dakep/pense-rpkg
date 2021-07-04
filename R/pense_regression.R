#' Compute (Adaptive) Elastic Net S-Estimates of Regression
#'
#' Compute elastic net S-estimates (PENSE estimates) along a grid of penalization levels with
#' optional penalty loadings for adaptive elastic net.
#'
#' @section Strategies for Using Starting Points:
#' The function supports several different strategies to compute, and use the provided starting
#' points for optimizing the PENSE objective function.
#'
#' Starting points are computed internally but can also be supplied via `other_starts`.
#' By default, starting points are computed internally by the EN-PY procedure for penalization
#' levels supplied in `enpy_lambda` (or the automatically generated grid of length `nlambda_enpy`).
#' By default, starting points computed by the EN-PY procedure are *shared* for all penalization
#' levels in `lambda` (or the automatically generated grid of length `nlambda`).
#' If the starting points should be *specific* to the penalization level the starting points'
#' penalization level, set the `enpy_specific` argument to `TRUE`.
#'
#' In addition to EN-PY initial estimates, the algorithm can also use the "0-based" strategy if
#' `add_zero_based = TRUE` (by default). Here, the 0-vector is used to start the optimization at
#' the largest penalization level in `lambda`. At subsequent penalization levels, the solution at
#' the previous penalization level is also used as starting point.
#'
#' At every penalization level, all starting points are explored using the loose numerical
#' tolerance `explore_tol`. Only the best `explore_solutions` are computed to the stringent
#' numerical tolerance `eps`.
#' Finally, only the best `max_solutions` are retained and carried forward as starting points for
#' the subsequent penalization level.
#'
#' @section Deprecated Arguments:
#' Starting with version 2.0.0, cross-validation is performed by separate function [pense_cv()].
#' Arguments related cross-validation cause an error when supplied to `pense()`.
#' Furthermore, the following arguments are deprecated as of version 2.0.0:
#' `initial`, `warm_reset`, `cl`, `options`, `init_options`, `en_options`.
#' If `pense()` is called with any of these arguments, warnings detail how to replace them.
#'
#' @param x `n` by `p` matrix of numeric predictors.
#' @param y vector of response values of length `n`.
#'          For binary classification, `y` should be a factor with 2 levels.
#' @param alpha elastic net penalty mixing parameter with \eqn{0 \le \alpha \le 1}.
#'    `alpha = 1` is the LASSO penalty, and `alpha = 0` the Ridge penalty.
#'    Can be a vector of several values, but `alpha = 0` cannot be mixed with other values.
#' @param nlambda number of penalization levels.
#' @param lambda_min_ratio Smallest value of the penalization level as a fraction of the largest
#'    level (i.e., the smallest value for which all coefficients are zero). The default depends on
#'    the sample size relative to the number of variables and `alpha`. If more observations than
#'    variables are available, the default is `1e-3 * alpha`, otherwise `1e-2 * alpha`.
#' @param nlambda_enpy number of penalization levels where the EN-PY initial estimate is computed.
#' @param penalty_loadings a vector of positive penalty loadings (a.k.a. weights) for different
#'    penalization of each coefficient. Only allowed for `alpha` > 0.
#' @param lambda optional user-supplied sequence of penalization levels. If given and not `NULL`,
#'    `nlambda` and `lambda_min_ratio` are ignored.
#' @param enpy_lambda optional user-supplied sequence of penalization levels at which EN-PY
#'    initial estimates are computed. If given and not `NULL`, `nlambda_enpy` is ignored.
#' @param other_starts a list of other staring points, created by [starting_point()].
#'    If the output of [enpy_initial_estimates()] is given, the starting points will be *shared*
#'    among all penalization levels.
#'    Note that if a the starting point is *specific* to a penalization level, this penalization
#'    level is added to the grid of penalization levels (either the manually specified grid in
#'    `lambda` or the automatically generated grid of size `nlambda`).
#'    If `standardize = TRUE`, the starting points are also scaled.
#' @param standardize logical flag to standardize the `x` variables prior to fitting the PENSE
#'    estimates. Coefficients are always returned on the original scale. This can fail for
#'    variables with a large proportion of a single value (e.g., zero-inflated data).
#'    In this case, either compute with `standardize = FALSE` or standardize the data manually.
#' @param intercept include an intercept in the model.
#' @param bdp desired breakdown point of the estimator, between 0 and 0.5. The actual
#'    breakdown point may be slightly larger/smaller to avoid instabilities of the S-loss.
#' @param cc tuning constant for the S-estimator. Default is to chosen based on the breakdown
#'   point \code{bdp}. Does *not* affect the estimated coefficients, only the estimated scale of
#'   the residuals.
#' @param eps numerical tolerance.
#' @param explore_solutions number of solutions to compute up to the desired precision `eps`.
#' @param explore_tol,explore_it numerical tolerance and maximum number of iterations for
#'    exploring possible solutions. The tolerance should be (much) looser than `eps` to be useful,
#'    and the number of iterations should also be much smaller than the maximum number of
#'    iterations given via `algorithm_opts`.
#' @param max_solutions only retain up to `max_solutions` unique solutions per penalization level.
#' @param comparison_tol numeric tolerance to determine if two solutions are equal.
#'    The comparison is first done on the absolute difference in the value of the objective
#'    function at the solution If this is less than `comparison_tol`, two solutions are deemed
#'    equal if the squared difference of the intercepts is less than `comparison_tol` and the
#'    squared \eqn{L_2} norm of the difference vector is less than `comparison_tol`.
#' @param add_zero_based also consider the 0-based regularization path. See details for a
#'    description.
#' @param enpy_specific use the EN-PY initial estimates only at the penalization level they
#'    are computed for. See details for a description.
#' @param sparse use sparse coefficient vectors.
#' @param ncores number of CPU cores to use in parallel. By default, only one CPU core is used.
#'    Not supported on all platforms, in which case a warning is given.
#' @param algorithm_opts options for the MM algorithm to compute the estimates.
#'    See [mm_algorithm_options()] for details.
#' @param mscale_opts options for the M-scale estimation. See [mscale_algorithm_options()]
#'    for details.
#' @param enpy_opts options for the ENPY initial estimates, created with the
#'    [enpy_options()] function. See [enpy_initial_estimates()] for details.
#' @param cv_k,cv_objective deprecated and ignored. See [pense_cv()] for estimating
#'    prediction performance via cross-validation.
#' @param ... ignored. See the section on deprecated parameters below.
#'
#' @return a list-like object with the following items
#'    \describe{
#'      \item{`alpha`}{the sequence of `alpha` parameters.}
#'      \item{`lambda`}{a list of sequences of penalization levels, one per `alpha` parameter.}
#'      \item{`estimates`}{a list of estimates. Each estimate contains the following information:
#'        \describe{
#'          \item{`intercept`}{intercept estimate.}
#'          \item{`beta`}{beta (slope) estimate.}
#'          \item{`lambda`}{penalization level at which the estimate is computed.}
#'          \item{`alpha`}{*alpha* hyper-parameter at which the estimate is computed.}
#'          \item{`bdp`}{chosen breakdown-point.}
#'          \item{`objf_value`}{value of the objective function at the solution.}
#'          \item{`statuscode`}{if `> 0` the algorithm experienced issues when
#'                              computing the estimate.}
#'          \item{`status`}{optional status message from the algorithm.}
#'        }
#'      }
#'      \item{`bdp`}{the actual breakdown point used.}
#'      \item{`call`}{the original call.}
#'    }
#'
#' @family functions to compute robust estimates
#' @seealso [pense_cv()] for selecting hyper-parameters via cross-validation.
#' @seealso [coef.pense_fit()] for extracting coefficient estimates.
#' @seealso [plot.pense_fit()] for plotting the regularization path.
#'
#' @example examples/pense_fit.R
#' @export
#' @aliases adapense
#' @importFrom lifecycle deprecated is_present deprecate_stop
pense <- function(x, y, alpha, nlambda = 50, nlambda_enpy = 10, lambda, lambda_min_ratio,
                  enpy_lambda, penalty_loadings, intercept = TRUE, bdp = 0.25, cc,
                  add_zero_based = TRUE, enpy_specific = FALSE, other_starts,
                  eps = 1e-6, explore_solutions = 10, explore_tol = 0.1, explore_it = 20,
                  max_solutions = 10, comparison_tol = sqrt(eps), sparse = FALSE,
                  ncores = 1, standardize = TRUE,
                  algorithm_opts = mm_algorithm_options(),
                  mscale_opts = mscale_algorithm_options(),
                  enpy_opts = enpy_options(), cv_k = deprecated(),
                  cv_objective = deprecated(), ...) {

  # Stop for CV-related options. Must migrate to `pense_cv`
  if (is_present(cv_k)) {
    deprecate_stop('2.0.0', 'pense(cv_k=)', 'pense_cv()')
  }
  if (is_present(cv_objective)) {
    deprecate_stop('2.0.0', 'pense(cv_objective=)', 'pense_cv()')
  }

  call <- match.call(expand.dots = TRUE)
  args <- as.list(call[-1L])
  args$standardize <- isTRUE(standardize)  # Ignore standardize = 'cv_only'!
  args <- do.call(.pense_args, args, envir = parent.frame())

  # Update BDP for numerical stability
  stable_bdp <- .find_stable_bdb_bisquare(
    n = length(args$std_data$y),
    desired_bdp = args$pense_opts$mscale$delta)
  args$pense_opts$mscale$delta <- stable_bdp

  # Call internal function
  fits <- mapply(
    args$alpha, args$lambda, args$enpy_lambda_inds,
    SIMPLIFY = FALSE, USE.NAMES = FALSE,
    FUN = function (alpha, lambda, enpy_lambda_inds) {
      fit <- .pense_internal(x = args$std_data$x, y = args$std_data$y,
                             alpha = alpha,
                             lambda = lambda,
                             enpy_lambda_inds = enpy_lambda_inds,
                             penalty_loadings = args$penalty_loadings,
                             pense_opts = args$pense_opts,
                             enpy_opts = args$enpy_opts,
                             optional_args = args$optional_args)

      # Retain only the best solution:
      fit$estimates <- lapply(fit$estimates, function (ests) {
        args$restore_coef_length(args$std_data$unstandardize_coefs(ests[[1L]]))
      })
      # Handle metrics
      fit$estimates <- .metrics_attrib(fit$estimates, fit$metrics)
      fit$lambda <- unlist(vapply(fit$estimates, FUN = `[[`, FUN.VALUE = numeric(1),
                                  'lambda'), use.names = FALSE, recursive = FALSE)
      fit$alpha <- alpha

      fit
    })

  structure(list(
    call = call,
    bdp = stable_bdp,
    lambda = lapply(fits, `[[`, 'lambda'),
    estimates = unlist(lapply(fits, `[[`, 'estimates'), recursive = FALSE, use.names = FALSE),
    alpha = vapply(fits, FUN.VALUE = numeric(1L), FUN = `[[`, 'alpha', USE.NAMES = FALSE)),
    class = c('pense', 'pense_fit'))
}

#' Cross-validation for (Adaptive) PENSE Estimates
#'
#' Perform (repeated) K-fold cross-validation for [pense()].
#'
#' @inheritParams pense
#' @param standardize whether to standardize the `x` variables prior to fitting the PENSE estimates.
#'    Can also be set to `"cv_only"`, in which case the input data is not standardized, but the
#'    training data in the CV folds is scaled to match the scaling of the input data.
#'    Coefficients are always returned on the original scale.
#'    This can fail for variables with a large proportion of a single value
#'    (e.g., zero-inflated data).
#'    In this case, either compute with `standardize = FALSE` or standardize the data manually.
#' @template cv_params
#' @inheritDotParams pense -standardize
#'
#' @seealso [pense()] for computing regularized S-estimates without cross-validation.
#'
#' @return a list-like object with the same components as returned by [pense()],
#'    plus the following:
#'    \describe{
#'      \item{`cvres`}{data frame of average cross-validated performance.}
#'    }
#'
#' @example examples/adapense_fit.R
#' @family functions to compute robust estimates with CV
#' @seealso [coef.pense_cvfit()] for extracting coefficient estimates.
#' @seealso [plot.pense_cvfit()] for plotting the CV performance or the regularization path.
#' @aliases adapense_cv
#' @export
#' @importFrom lifecycle deprecate_warn deprecated is_present
#' @importFrom stats sd
#' @importFrom rlang abort
pense_cv <- function(x, y, standardize = TRUE, lambda, cv_k, cv_repl = 1,
                     cv_metric = c('tau_size', 'mape', 'rmspe', 'auroc'), fit_all = TRUE,
                     cl = NULL, ...) {
  call <- match.call(expand.dots = TRUE)
  args <- do.call(.pense_args, as.list(call[-1L]), envir = parent.frame())

  fit_ses <- if (is.character(fit_all)) {
    unique(vapply(fit_all, FUN = .parse_se_string, FUN.VALUE = numeric(1L), only_fact = TRUE,
                  USE.NAMES = FALSE))
  } else if (isFALSE(fit_all)) {
    .parse_se_string('min', only_fact = TRUE)
  } else {
    TRUE
  }

  if (length(fit_ses) < 1L) {
    fit_ses <- TRUE
  }

  cv_k <- .as(cv_k[[1L]], 'integer')
  cv_repl <- .as(cv_repl[[1L]], 'integer')

  if (cv_k < 2L) {
    abort("`cv_k` must be greater than 1.")
  }

  if (cv_repl < 1L) {
    abort("`cv_repl` must be greater than 0.")
  }

  if (identical(cv_repl, 1L) && any(fit_ses > 0)) {
    warn("To use `fit_all = \"se\"`, `cv_repl` must be 2 or greater.")
    fit_ses <- 0
  }

  if (!is.null(cl)) {
    if (args$pense_opts$num_threads > 1L) {
      abort("`cl` can only be used if `ncores = 1`.")
    }
    if (!is(cl, 'cluster')) {
      abort("`cl` must be a valid `parallel` cluster.")
    }
  }

  cv_metric <- if (is.null(call$cv_metric) && args$binary_response) {
    cv_measure_str <- 'auroc'
    .cv_auroc
  } else if (is.character(cv_metric)) {
    cv_measure_str <- match.arg(cv_metric)
    switch(cv_measure_str, mape = .cv_mape, rmspe = .cv_rmspe, tau_size = tau_size,
           auroc = if (args$binary_response) {
             .cv_auroc
           } else {
             abort("cv_metric=\"auroc\" is only valid for binary responses.")
           })
  } else {
    cv_measure_str <- 'user_fun'
    match.fun(cv_metric)
  }

  if (is.null(formals(cv_metric))) {
    abort("Function `cv_metric` must accept at least 1 argument.")
  }

  # Get a common seed to be used for every alpha value
  fit_seed <- sample.int(.Machine$integer.max, 1L)

  cv_curves <- mapply(
    args$alpha, args$lambda, args$enpy_lambda_inds,
    SIMPLIFY = FALSE, USE.NAMES = FALSE,
    FUN = function (alpha, lambda, enpy_lambda_inds) {
      cv_fun <- function (train_data, test_ind) {
        stable_bdp <- .find_stable_bdb_bisquare(
          n = length(train_data$y),
          desired_bdp = args$pense_opts$mscale$delta)
        args$pense_opts$mscale$delta <- stable_bdp

        cv_fit <- .pense_internal(x = train_data$x, y = train_data$y,
                                  alpha = alpha,
                                  lambda = lambda,
                                  enpy_lambda_inds = enpy_lambda_inds,
                                  penalty_loadings = args$penalty_loadings,
                                  pense_opts = args$pense_opts,
                                  enpy_opts = args$enpy_opts,
                                  optional_args = args$optional_args)

        # Return only best local optima
        lapply(cv_fit$estimates, `[[`, 1L)
      }

      set.seed(fit_seed)
      cv_perf <- .run_replicated_cv(args$std_data,
                                    cv_k = cv_k,
                                    cv_repl = cv_repl,
                                    metric = cv_metric,
                                    cv_est_fun = cv_fun,
                                    par_cluster = cl)
      data.frame(lambda = lambda, alpha = alpha,
                 cvavg = rowMeans(cv_perf),
                 cvse = if (cv_repl > 1L) { apply(cv_perf, 1, sd) } else { 0 })
    })

  cv_curves <- do.call(rbind, cv_curves)

  if (isTRUE(fit_ses)) {
    fit_lambda <- args$lambda
    fit_enpy_lambda_inds <- args$enpy_lambda_inds
  } else {
    fit_lambda <- lapply(args$alpha, function (alpha) {
      rows <- which((cv_curves$alpha - alpha)^2 < .Machine$double.eps)

      lambda_inds <- vapply(fit_ses, FUN.VALUE = numeric(1L), FUN = function (se_fact) {
        which(.cv_se_selection(cv_curves$cvavg[rows], cv_curves$cvse[rows], se_fact) == 'se_fact')
      })
      unique(cv_curves$lambda[rows[lambda_inds]])
    })
    fit_enpy_lambda_inds <- lapply(fit_lambda, seq_along)
  }

  stable_bdp <- .find_stable_bdb_bisquare(
    n = length(args$std_data$y),
    desired_bdp = args$pense_opts$mscale$delta)
  args$pense_opts$mscale$delta <- stable_bdp

  # Call internal function
  fits <- mapply(
    args$alpha, fit_lambda, fit_enpy_lambda_inds,
    SIMPLIFY = FALSE, USE.NAMES = FALSE,
    FUN = function (alpha, lambda, enpy_lambda_inds) {
      fit <- .pense_internal(x = args$std_data$x, y = args$std_data$y,
                             alpha = alpha,
                             lambda = lambda,
                             enpy_lambda_inds = enpy_lambda_inds,
                             penalty_loadings = args$penalty_loadings,
                             pense_opts = args$pense_opts,
                             enpy_opts = args$enpy_opts,
                             optional_args = args$optional_args)

      # Retain only the best solution:
      fit$estimates <- lapply(fit$estimates, function (ests) {
        args$restore_coef_length(args$std_data$unstandardize_coefs(ests[[1L]]))
      })
      # Handle metrics
      fit$estimates <- .metrics_attrib(fit$estimates, fit$metrics)
      fit$lambda <- unlist(vapply(fit$estimates, FUN = `[[`, FUN.VALUE = numeric(1),
                                  'lambda'), use.names = FALSE, recursive = FALSE)
      fit$alpha <- alpha
      fit
    })

  structure(list(
    call = call,
    bdp = stable_bdp,
    lambda = lapply(fits, `[[`, 'lambda'),
    alpha = vapply(fits, FUN.VALUE = numeric(1L), FUN = `[[`, 'alpha', USE.NAMES = FALSE),
    cvres = cv_curves,
    cv_measure = cv_measure_str,
    estimates = unlist(lapply(fits, `[[`, 'estimates'), recursive = FALSE, use.names = FALSE)),
    class = c('pense', 'pense_cvfit'))
}

#' @description `adapense_cv()` is a convenience wrapper to compute adaptive PENSE estimates.
#'
#' @details
#' `adapense_cv()` is a convenience wrapper which performs 3 steps:
#'
#' 1. compute preliminary estimates via `pense_cv(..., alpha = alpha_preliminary)`,
#' 2. computes the penalty loadings from the estimate `beta` with best prediction performance by
#'    `adapense_loadings = 1 / abs(beta)^exponent`, and
#' 3. compute the adaptive PENSE estimates via `pense_cv(..., penalty_loadings = adapense_loadings)`.
#'
#' @param alpha_preliminary `alpha` parameter for the preliminary estimate.
#' @param exponent the exponent for computing the penalty loadings based on the preliminary estimate.
#'
#' @return a list-like object as returned by [pense_cv()] plus the following
#'    \describe{
#'      \item{`preliminary`}{the CV results for the preliminary estimate.}
#'      \item{`exponent`}{exponent used to compute the penalty loadings.}
#'      \item{`penalty_loadings`}{the penalty loadings used for the adaptive PENSE estimate.}
#'    }
#'
#' @family functions to compute robust estimates with CV
#' @rdname pense_cv
#' @export
#' @importFrom stats coef
adapense_cv <- function (x, y, alpha, alpha_preliminary = 0, exponent = 1, ...) {
  call <- match.call(expand.dots = TRUE)
  if (!is.null(call$penalty_loadings)) {
    abort(paste("Argument `penalty_loadings` not valid for `adapense_cv`.",
                "Penalty loadings are determined internally."))
  }
  exponent <- .as(exponent[[1L]], 'numeric')

  # Compute preliminary estimate
  prelim_call <- call
  prelim_call[[1L]] <- quote(pense::pense_cv)
  prelim_call$alpha <- .as(alpha_preliminary[[1]], 'numeric')
  prelim_call$alpha_preliminary <- NULL
  prelim_call$exponent <- NULL
  prelim <- eval.parent(prelim_call)

  prelim_coef <- coef(prelim, sparse = FALSE, concat = FALSE)
  pen_loadings <- abs(prelim_coef$beta)^(-exponent)

  adapense <- pense_cv(x, y, alpha = alpha, penalty_loadings = pen_loadings, ...)
  adapense$call <- call
  adapense$exponent <- exponent
  adapense$preliminary <- prelim
  adapense$penalty_loadings <- pen_loadings
  class(adapense) <- c('adapense', class(adapense))
  return(adapense)
}

## Make a list of initial estimates
##
## @return a list with 2 components:
##   `extended_lambda` an extended grid of penalization levels to contain both the given
##                     lambda values plus the lambda values in `other_starts`
##   `starting_points` a list the same length as `extended_lambda` with a list of initial estimates
##                     for each value in `extended_lambda`.
#' @importFrom rlang warn
.make_initest_list <- function (other_starts, lambda, alpha, sparse) {
  if (length(other_starts) == 0L) {
    return(list(extended_lambda = lambda, starting_points = rep.int(list(list()), length(lambda))))
  }

  # Check for wrong starting points without lambda
  init_est_lambda <- unlist(lapply(other_starts, `[[`, 'lambda'),
                            use.names = FALSE, recursive = FALSE)

  if (length(init_est_lambda) != length(other_starts)) {
    abort(paste("Some starting points in `other_starts` are marked as \"specific\"",
                "but do not have a `lambda` component."))
  }

  init_est_lambda <- .as(init_est_lambda, 'numeric')

  # Check for wrong starting points without alpha
  init_est_alpha <- unlist(lapply(other_starts, `[[`, 'alpha'),
                           use.names = FALSE, recursive = FALSE)

  if (length(init_est_alpha) != length(other_starts)) {
    abort(paste("Some starting points in `other_starts` are marked as \"specific\"",
                "but do not have an `alpha` component."))
  }

  init_est_alpha <- .as(init_est_alpha, 'numeric')
  correct_alpha <- which(abs(init_est_alpha - alpha) < sqrt(.Machine$double.eps))

  if (length(correct_alpha) == 0L) {
    return(list(extended_lambda = lambda, starting_points = rep.int(list(list()), length(lambda))))
  }
  init_est_lambda <- init_est_lambda[correct_alpha]
  other_starts <- other_starts[correct_alpha]

  init_est_inds <- .approx_match(init_est_lambda, lambda)
  new_initest_lambda <- which(is.na(init_est_inds))

  if (length(new_initest_lambda) > 0L) {
    # Some starts are for unknown lambda. Add lambdas to the grid!
    lambda <- sort(c(lambda, unique(init_est_lambda[new_initest_lambda])), decreasing = TRUE)
    init_est_inds <- .approx_match(init_est_lambda, lambda)
  }

  starting_points <- lapply(seq_along(lambda), function (i) {
    matches <- which(i == init_est_inds)
    if (length(matches) > 0L) {
      return(other_starts[matches])
    } else {
      return(list())
    }
  })

  return(list(extended_lambda = lambda, starting_points = starting_points))
}

## Get the smallest lambda such that the PENSE estimate gives the empty model.
.pense_max_lambda <- function (x, y, alpha, pense_options, penalty_loadings = NULL) {
  optional_args <- list()
  if (!is.null(penalty_loadings)) {
    optional_args$pen_loadings <- penalty_loadings
  }
  .Call(C_pense_max_lambda, x, y, pense_options, optional_args) / max(0.01, alpha)
}

## Generate a log-spaced grid of decreasing lambda values
.pense_lambda_grid <- function (x, y, alpha, nlambda, lambda_min_ratio, pense_options,
                                penalty_loadings) {
  alpha <- max(0.01, alpha)
  x_dim <- dim(x)
  if (is.null(lambda_min_ratio)) {
    lambda_min_ratio <- alpha * if (x_dim[[1L]] > x_dim[[2L]]) { 1e-3 } else { 1e-2 }
  }
  max_lambda <- .pense_max_lambda(x, y, alpha, pense_options, penalty_loadings)
  rev(exp(seq(log(lambda_min_ratio * max_lambda), log(max_lambda), length.out = nlambda)))
}

## Perform some final input adjustments and call the internal C++ code.
.pense_internal <- function(x, y, alpha, lambda, enpy_lambda_inds, penalty_loadings = NULL,
                            pense_opts, enpy_opts, optional_args) {
  # Create penalties-list, without sorting the lambda sequence
  penalties <- lapply(lambda, function (l) { list(lambda = l, alpha = alpha) })

  if (!is.null(penalty_loadings)) {
    optional_args$pen_loadings <- penalty_loadings
  }

  .Call(C_pense_regression, x, y, penalties, enpy_lambda_inds, pense_opts, enpy_opts, optional_args)
}

#' @importFrom lifecycle deprecated is_present deprecate_stop
#' @importFrom rlang warn abort
#' @importFrom methods is
#' @importFrom stats runif
.pense_args <- function (x, y, alpha, nlambda = 50, nlambda_enpy = 10, lambda, lambda_min_ratio,
                         enpy_lambda, penalty_loadings, intercept = TRUE, bdp = 0.25, cc = NULL,
                         add_zero_based = TRUE, enpy_specific = FALSE, other_starts,
                         eps = 1e-6, explore_solutions = 10, explore_tol = 0.1, explore_it = 20,
                         max_solutions = 10, comparison_tol = sqrt(eps), sparse = FALSE,
                         ncores = 1, standardize = TRUE,
                         algorithm_opts = mm_algorithm_options(),
                         mscale_opts = mscale_algorithm_options(),
                         enpy_opts = enpy_options(),
                         options = deprecated(), init_options = deprecated(),
                         en_options = deprecated(),
                         initial = deprecated(), warm_reset = deprecated(), ...) {
  args_call <- match.call(expand.dots = FALSE)
  optional_args <- list()

  ## Check presence of deprecated arguments
  # Translate options for initial estimates
  if (is_present(warm_reset)) {
    if (is.null(args_call$nlambda_enpy)) {
      deprecate_warn('2.0.0', 'pense(warm_reset=)', 'pense(nlambda_enpy=)')
      nlambda_enpy <- warm_reset
    } else {
      deprecate_warn('2.0.0', 'pense(warm_reset=)',
                     details = paste("Superseding argument `nlambda_enpy` is also present",
                                     "and will be used instead of `warm_reset`."))
    }
  }
  if (is_present(initial)) {
    if (!is.null(args_call$nlambda_enpy)) {
      abort(paste("The `initial` argument of `pense()` is deprecated as of pense 2.0.0",
                  "and conflicts with the provided `nlambda_enpy` argument."))
    }
    if (!is.null(args_call$enpy_lambda)) {
      abort(paste("The `initial` argument of `pense()` is deprecated as of pense 2.0.0",
                  "and conflicts with the provided `enpy_lambda` argument."))
    }

    deprecate_warn('2.0.0', 'pense(initial=)',
                   details = paste("Please use arguments `nlambda_enpy`, `enpy_specific`,",
                                   "`add_zero_based` and `other_starts` instead."))
    initial <- match.arg(initial, c('warm', 'cold'))
    add_zero_based <- TRUE
    enpy_specific <- TRUE
    if (initial == 'cold') {
      nlambda_enpy <- if (missing(lambda)) { nlambda } else { length(lambda) }
    }
  }

  # Translate algorithm options
  if (is_present(en_options)) {
    deprecate_warn('2.0.0', 'pense(en_options=)',
                   details = paste("Please specify the LS-EN algorithm with arguments",
                                   "`algorithm_opts` and `enpy_opts`."))
    if (is.null(args_call$algorithm_opts)) {
      algorithm_opts$en_opts <- en_options
    }
  }
  if (is_present(init_options)) {
    if (is.null(args_call$enpy_opts)) {
      deprecate_warn('2.0.0', 'pense(init_options=)', 'pense(enpy_opts=)')
      enpy_opts <- enpy_options(max_it = init_options$maxit,
                                keep_psc_proportion = init_options$keepPSCProportion,
                                keep_residuals_measure = init_options$keepResidualsMethod,
                                keep_residuals_proportion = init_options$keepResidualsProportion,
                                keep_residuals_threshold = init_options$keepResidualsThreshold,
                                en_algorithm_opts = if (is_present(en_options)) {
                                  en_options
                                } else {
                                  NULL
                                })
      explore_solutions <- init_options$keepSolutions
    } else {
      deprecate_warn('2.0.0', 'pense(init_options=)',
                     details = paste("Superseding argument `enpy_opts` is also present and will",
                                     "be used instead of `init_options`."))
    }
  }
  if (is_present(options)) {
    deprecate_warn('2.0.0', 'pense(options=)',
                   details = paste("Please use arguments `bdp`, `eps`, `mscale_opts`",
                                   "and `algorithm_opts` to `pense()` instead"))
    if (is.null(args_call$bdp)) { bdp <- options$bdp }
    if (is.null(args_call$cc)) { cc <- options$cc }
    if (is.null(args_call$eps)) { eps <- options$eps }
    if (is.null(args_call$mscale_opts)) { mscale_opts <- options$mscale_opts }
    if (is.null(args_call$algorithm_opts) && !is.null(algorithm_opts$max_it)) {
      algorithm_opts$max_it <- options$maxit
    }
  }

  ## Process input arguments
  response <- .validate_response(y)
  y <- response$values
  x_dim <- dim(x)

  if (length(y) != x_dim[[1L]]) {
    abort("Number of observations in `x` and `y` does not match.")
  } else if (x_dim[[2L]] <= 1L) {
    abort("`x` must be a matrix with at least 2 columns.")
  }

  alpha <- .as(alpha, 'numeric')
  if (any(alpha < 0 | alpha > 1)) {
    abort("`alpha` is outside 0 and 1.")
  } else if (any(alpha < sqrt(.Machine$double.eps))) {
    alpha[which(alpha < sqrt(.Machine$double.eps))] <- 0
    if (any(alpha > 0)) {
      abort("`alpha=0` cannot be mixed with other `alpha` values.")
    }
  }

  if (!missing(enpy_lambda)) {
    nlambda_enpy <- length(enpy_lambda)
  }

  pense_opts <- list(algo_opts = algorithm_opts,
                     strategy_0 = isTRUE(add_zero_based),
                     strategy_enpy_individual = isTRUE(enpy_specific) && (nlambda_enpy > 0L),
                     strategy_enpy_shared = !isTRUE(enpy_specific) && (nlambda_enpy > 0L),
                     strategy_other_individual = FALSE,
                     strategy_other_shared = FALSE,
                     algorithm = .pense_algorithm_id(algorithm_opts),
                     intercept = !isFALSE(intercept),
                     eps = .as(eps[[1L]], 'numeric'),
                     comparison_tol = .as(comparison_tol[[1L]], 'numeric'),
                     explore_tol = .as(explore_tol[[1L]], 'numeric'),
                     explore_it = .as(explore_it[[1L]], 'integer'),
                     nr_tracks = .as(explore_solutions[[1L]], 'integer'),
                     max_optima = .as(max_solutions[[1L]], 'integer'),
                     num_threads = max(1L, .as(ncores[[1L]], 'integer')),
                     sparse = isTRUE(sparse),
                     mscale = .full_mscale_algo_options(bdp = bdp, cc = cc,
                                                        mscale_opts = mscale_opts))

  if (pense_opts$explore_tol < pense_opts$eps) {
    abort("`explore_tol` must not be less than `eps`")
  }
  if (pense_opts$comparison_tol < pense_opts$eps) {
    abort("`comparison_tol` must not be less than `eps`")
  }
  if (pense_opts$explore_it < 1L) {
    abort("`explore_it` must not be less than 0")
  }

  # Check EN algorithm for ENPY
  enpy_opts$en_options <- .select_en_algorithm(enpy_opts$en_options, alpha, pense_opts$sparse, eps)
  pense_opts$sparse <- enpy_opts$en_options$sparse

  # If using the MM algorithm, ensure that the EN options are set.
  if (pense_opts$algorithm == 1L) {
    pense_opts$algo_opts$en_options <- .select_en_algorithm(pense_opts$algo_opts$en_options, alpha,
                                                            pense_opts$sparse, eps)
    if (!isTRUE(pense_opts$sparse == pense_opts$algo_opts$en_options$sparse)) {
      abort("The `sparse` option for the EN-PY algorithm and the MM algorithm for PENSE disagree.")
    }
  }

  # Set the number of cores for the ENPY options
  if (pense_opts$num_threads > 1L && !isTRUE(.k_multithreading_support)) {
    warn("Multithreading not supported. Using only 1 core.")
    pense_opts$num_threads <- 1L
  }
  enpy_opts$num_threads <- pense_opts$num_threads

  # Standardizing the data
  standardize <- if (is.character(standardize)) {
    if (pmatch(standardize[[1L]], 'cv_only', nomatch = 0L) == 1L) {
      standardize <- 'cv_only'
    } else {
      abort("`standardize` must be either TRUE/FALSE or \"cv_only\".")
    }
  } else {
    isTRUE(standardize)
  }

  # Check penalty loadings
  if (!missing(penalty_loadings) && !is.null(penalty_loadings)) {
    checked_pls <- .prepare_penalty_loadings(penalty_loadings, x = x, alpha = alpha,
                                             sparse = pense_opts$sparse)
    penalty_loadings <- checked_pls$loadings
    restore_coef_length <- checked_pls$restore_fun
    x <- checked_pls$trimmed_x
  } else {
    restore_coef_length <- function (coef) coef
    penalty_loadings <- NULL
  }

  if (ncol(x) == 0L) {
    pense_opts$intercept <- TRUE
    warn("All values in `penalty_loadings` are infinite. Only computing the intercept.")
    std_data <- .standardize_data(matrix(runif(x_dim[[1L]]), ncol = 1L), y, intercept = TRUE,
                                  sparse = pense_opts$sparse, standardize = standardize,
                                  robust = TRUE, mscale_opts = mscale_opts,
                                  bdp = pense_opts$mscale$delta, scale_cc = pense_opts$mscale$cc)
    # Compute only the 0-based solution.
    pense_opts$strategy_enpy_individual <- FALSE
    pense_opts$strategy_enpy_shared <- FALSE
    pense_opts$strategy_0 <- TRUE
    lambda <- lapply(alpha, FUN = .pense_lambda_grid,
                     x = std_data$x, y = std_data$y, nlambda = 1, lambda_min_ratio = 1,
                     pense_options = pense_opts, penalty_loadings = NULL)

    return(list(std_data = std_data,
                alpha = alpha,
                lambda = lambda,
                enpy_lambda_inds = integer(0L),
                penalty_loadings = NULL,
                pense_opts = pense_opts,
                enpy_opts = enpy_opts,
                optional_args = optional_args,
                restore_coef_length = restore_coef_length))
  }

  std_data <- .standardize_data(x, y, intercept = pense_opts$intercept, standardize = standardize,
                                robust = TRUE, sparse = pense_opts$sparse,
                                mscale_opts = mscale_opts, bdp = pense_opts$mscale$delta,
                                scale_cc = pense_opts$mscale$cc)

  # Scale penalty loadings appropriately
  penalty_loadings <- penalty_loadings / std_data$scale_x
  if (length(penalty_loadings) == 0L) {
    penalty_loadings <- NULL
  }

  # Determine lambda grid
  lambda <- if (missing(lambda) || is.null(lambda)) {
    if (missing(lambda_min_ratio)) {
      lambda_min_ratio <- NULL
    }
    lapply(alpha, FUN = .pense_lambda_grid,
           x = std_data$x, y = std_data$y, nlambda = nlambda,
           lambda_min_ratio = lambda_min_ratio, pense_options = pense_opts,
           penalty_loadings = penalty_loadings)
  } else if (!is.list(lambda)) {
    rep.int(list(sort(.as(lambda, 'numeric'), decreasing = TRUE)), length(alpha))
  } else if (identical(length(lambda), length(alpha))) {
    lapply(lambda, function (l) {
      sort(.as(l, 'numeric'), decreasing = TRUE)
    })
  } else {
    abort("`lambda` must either be a numeric vector or a list the same length as `alpha`.")
  }

  # Split the `other_starts` into individual and shared starts.
  if (!missing(other_starts)) {
    if (is(other_starts, 'starting_point')) {
      other_starts <- structure(list(other_starts), class = 'starting_points')
    } else if (!is(other_starts, 'starting_points')) {
      abort(paste("`other_starts` must be a list of starting points created by",
                  "`starting_point()`, `enpy_initial_estimates()`, or a combination thereof."))
    }

    # Identify which other starts are shared and which are specific.
    other_starts_shared <- vapply(other_starts, FUN.VALUE = logical(1L),
                                  FUN = is, 'shared_starting_point')
    other_starts_specific <- vapply(other_starts, FUN.VALUE = logical(1L),
                                    FUN = is, 'specific_starting_point')

    # Ensure the `beta` coefficients in `other_starts` agree with the desired vector class
    # (sparse vs. dense) and standardize them.
    other_starts <- lapply(.sparsify_other_starts(other_starts, pense_opts$sparse),
                           std_data$standardize_coefs)

    if (any(other_starts_shared)) {
      pense_opts$strategy_other_shared <- TRUE
      optional_args$shared_starts <- other_starts[other_starts_shared]
    }
    if (any(other_starts_specific)) {
      pense_opts$strategy_other_individual <- TRUE
      ind_starts <- list()
      for (ai in seq_along(alpha)) {
        new_ind_starts <- .make_initest_list(other_starts[other_starts_specific], lambda[[ai]],
                                             alpha = alpha[[ai]], sparse = pense_opts$sparse)
        ind_starts <- c(ind_starts, new_ind_starts$starting_points)
        lambda[[ai]] <- new_ind_starts$extended_lambda
      }
      optional_args$individual_starts <- ind_starts
    }
  }

  # Determine ENPY lambda grid
  enpy_lambda_inds <- if (pense_opts$strategy_enpy_individual || pense_opts$strategy_enpy_shared) {
    if (missing(enpy_lambda) || is.null(enpy_lambda)) {
      lapply(lambda, function (l) {
        nlambda_enpy <- min(length(l), nlambda_enpy)
        as.integer(ceiling(seq(1, length(l), length.out = nlambda_enpy + 1))[-(nlambda_enpy + 1)])
      })
    } else if (is.list(enpy_lambda)) {
      mapply(enpy_lambda, lambda, SIMPLIFY = FALSE, USE.NAMES = FALSE, FUN = .approx_match)
    } else if (is.numeric(enpy_lambda)) {
      vapply(lambda, FUN.VALUE = integer(length(enpy_lambda)),
             FUN = .approx_match, x = enpy_lambda)
    } else {
      abort("`enpy_lambda` must either be a numeric vector or a list the same length as `alpha`.")
    }
  } else {
    integer(0L)
  }

  # Extend lambda grids if necessary
  for (ai in seq_along(alpha)) {
    if (anyNA(enpy_lambda_inds[[ai]])) {
      if (is.list(enpy_lambda)) {
        lambda[[ai]] <- sort(c(lambda[[ai]], enpy_lambda[[ai]]), decreasing = TRUE)
        enpy_lambda_inds[[ai]] <- .approx_match(enpy_lambda[[ai]], lambda[[ai]])
      } else if (is.numeric(enpy_lambda)) {
        lambda[[ai]] <- sort(c(lambda[[ai]], enpy_lambda), decreasing = TRUE)
        enpy_lambda_inds[[ai]] <- .approx_match(enpy_lambda, lambda[[ai]])
      }
    }
  }

  for (lambda_grid in lambda) {
    if (any(lambda_grid < .Machine$double.eps)) {
      abort("All values in `lambda` must be positive.")
    }
  }

  return(list(std_data = std_data,
              binary_response = response$binary,
              alpha = alpha,
              lambda = lambda,
              enpy_lambda_inds = enpy_lambda_inds,
              penalty_loadings = penalty_loadings,
              pense_opts = pense_opts,
              enpy_opts = enpy_opts,
              optional_args = optional_args,
              restore_coef_length = restore_coef_length))
}
