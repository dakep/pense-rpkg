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
#' @param bdp desired breakdown point of the estimator, between 0.05 and 0.5. The actual
#'    breakdown point may be slightly larger/smaller to avoid instabilities of the S-loss.
#' @param cc tuning constant for the S-estimator. Default is chosen based on the breakdown
#'   point \code{bdp}. This affects the estimated coefficients only if
#'   `standardize=TRUE`. Otherwise only the estimated scale of the residuals
#'   would be affected.
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
#' @param carry_forward carry the best solutions forward to the next penalty
#'   level.
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
pense <- function(x, y, alpha, nlambda = 50, nlambda_enpy = 10, lambda,
                  lambda_min_ratio, enpy_lambda, penalty_loadings,
                  intercept = TRUE, bdp = 0.25, cc,
                  add_zero_based = TRUE, enpy_specific = FALSE, other_starts,
                  carry_forward = TRUE,
                  eps = 1e-6, explore_solutions = 10, explore_tol = 0.1,
                  explore_it = 5, max_solutions = 1,
                  comparison_tol = sqrt(eps), sparse = FALSE,
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
  call$standardize <- isTRUE(standardize)

  args_env <- new.env(parent = parent.frame())
  args_env$pense <- .pense_args
  call[[1]] <- quote(pense)
  args <- eval(call, envir = args_env)

  # Update BDP for numerical stability
  stable_bdp <- .find_stable_bdb_bisquare(
    n = length(args$std_data$y),
    desired_bdp = args$pense_opts$mscale$delta)
  args$pense_opts$mscale$delta <- stable_bdp

  # Call internal function
  fits <- .pense_internal_multi(args)

  structure(list(
    call = match.call(expand.dots = TRUE),
    bdp = stable_bdp,
    lambda = args$lambda,
    metrics = lapply(fits, function (f) { attr(f$estimates, 'metrics') }),
    estimates = unlist(lapply(fits, `[[`, 'estimates'), recursive = FALSE,
                       use.names = FALSE),
    alpha = vapply(fits, FUN.VALUE = numeric(1L), FUN = `[[`, 'alpha',
                   USE.NAMES = FALSE)),
    class = c('pense', 'pense_fit'))
}

#' Cross-validation for (Adaptive) PENSE Estimates
#'
#' Perform (repeated) K-fold cross-validation for [pense()].
#'
#' @inheritParams pense
#' @param standardize whether to standardize the `x` variables prior to fitting
#'    the PENSE estimates. Can also be set to `"cv_only"`, in which case the
#'    input data is not standardized, but the training data in the CV folds is
#'    scaled to match the scaling of the input data.
#'    Coefficients are always returned on the original scale.
#'    This can fail for variables with a large proportion of a single value
#'    (e.g., zero-inflated data).
#'    In this case, either compute with `standardize = FALSE` or standardize
#'    the data manually.
#' @param fold_starts how to determine starting values in the
#'    cross-validation folds. If `"full"` (default), use the best solution from
#'    the fit to the full data as starting value. This implies
#'    `fit_all=TRUE`.
#'    If `"enpy"` compute separate ENPY initial estimates in each fold.
#'    The option `"both"` uses both.
#'    These starts are in addition to the starts provided in `other_starts`.
#' @template cv_params
#' @inheritDotParams pense -standardize
#'
#' @seealso [pense()] for computing regularized S-estimates without
#'   cross-validation.
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
#' @seealso [plot.pense_cvfit()] for plotting the CV performance or the
#'    regularization path.
#' @aliases adapense_cv
#' @export
#' @importFrom lifecycle deprecate_stop deprecated is_present
#' @importFrom stats sd
#' @importFrom rlang abort
pense_cv <- function(x, y, standardize = TRUE, lambda, cv_k, cv_repl = 1,
                     cv_metric = c('tau_size', 'mape', 'rmspe', 'auroc'),
                     fit_all = TRUE,
                     fold_starts = c('full', 'enpy', 'both'),
                     cl = NULL, ...) {
  call <- match.call(expand.dots = TRUE)
  args_env <- new.env(parent = parent.frame())
  args_env$pense_cv <- .pense_args
  call[[1]] <- quote(pense_cv)
  args <- eval(call, envir = args_env)

  fit_ses <- if (is.character(fit_all)) {
    unique(vapply(fit_all, FUN = .parse_se_string, FUN.VALUE = numeric(1L),
                  only_fact = TRUE, USE.NAMES = FALSE))
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

  if (identical(cv_repl, 1L) && !isTRUE(fit_ses) && any(fit_ses > 0)) {
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

  n_total <- length(args$std_data$y)
  n_in_folds <- n_total - ceiling(n_total / cv_k)
  if (!isTRUE(args$pense_opts$mscale$delta * n_total / n_in_folds <= 0.5)) {
    abort(paste("The desired breakdown point cannot be achieved in",
                "all CV folds. Either increase the number of CV folds",
                "or decrease the desired breakdown point."))
  }

  fold_starts <- match.arg(fold_starts)

  if (identical(fold_starts, 'full') || identical(fold_starts, 'both')) {
    fit_ses <- TRUE
  }

  # Determine "stable" bdp
  stable_bdp <- .find_stable_bdb_bisquare(
    n = length(args$std_data$y),
    desired_bdp = args$pense_opts$mscale$delta)
  args$pense_opts$mscale$delta <- stable_bdp

  # If we need to estimate solutions for all lambda, do it prior to CV.
  fits <- list()
  other_starts <- NULL

  if (isTRUE(fit_ses)) {
    fits <- .pense_internal_multi(args)
  }

  other_starts <- if (isTRUE(fit_ses)) {
    make_other_starts <- function (fit_alpha, lambda_seq) {
      # If there are other individual starts, only use the ones with
      # correct `alpha`
      n_user_ind_starts <- length(args$optional_args$individual_starts)
      old_starts <- if (n_user_ind_starts > 0L) {
        .filter_list(args$optional_args$individual_starts, 'alpha',
                     fit_alpha$alpha)
      } else {
        list()
      }

      std_ests <- lapply(fit_alpha$estimates, function (est) {
        est$beta <- est$std_beta
        est$intercept <- est$std_intercept
        est$std_beta <- NULL
        est$std_intercept <- NULL
        est
      })

      .make_initest_list(c(old_starts, std_ests),
                         lambda = lambda_seq,
                         alpha = fit_alpha$alpha,
                         sparse = args$pense_opts$sparse)$starting_points
    }

    mapply(fits, args$lambda, SIMPLIFY = FALSE,
           FUN = make_other_starts)
  } else {
    lapply(args$alpha, function (alpha) {
      if (length(args$optional_args$individual_starts) > 0L) {
        .filter_list(args$optional_args$individual_starts, 'alpha', alpha)
      } else {
        list()
      }
    })
  }

  # Get a common seed to be used for every alpha value
  fit_seed <- sample.int(.Machine$integer.max, 1L)

  cv_est_fun <- function (train_data, test_ind, handler_args) {
    # Determine stable bdp separately for this fold
    desired_bdp <- handler_args$args$pense_opts$mscale$delta *
      (length(train_data$y) + length(test_ind)) / length(train_data$y)

    stable_bdp <- .find_stable_bdb_bisquare(
      n = length(train_data$y),
      desired_bdp = desired_bdp)

    handler_args$args$pense_opts$mscale$delta <- stable_bdp

    cv_fit <- .pense_internal(
      x = train_data$x,
      y = train_data$y,
      alpha = handler_args$alpha,
      lambda = handler_args$lambda,
      enpy_lambda_inds = handler_args$enpy_lambda_inds,
      penalty_loadings = handler_args$args$penalty_loadings,
      pense_opts = handler_args$args$pense_opts,
      enpy_opts = handler_args$args$enpy_opts,
      optional_args = handler_args$args$optional_args)

    # Return only best local optima
    lapply(cv_fit$estimates, `[[`, 1L)
  }

  cv_est_dispatch <- function (alpha, lambda, enpy_lambda_inds,
                               other_starts) {
    handler_args <- list(alpha = alpha,
                         lambda = lambda,
                         cv_k = cv_k,
                         enpy_lambda_inds = enpy_lambda_inds,
                         args = args)

    if (identical(fold_starts, 'full')) {
      handler_args$enpy_lambda_inds <- integer(0L)
    }

    if (!identical(fold_starts, 'enpy')) {
      handler_args$args$pense_opts$strategy_other_individual <- TRUE
    }

    handler_args$args$optional_args$individual_starts <- other_starts

    set.seed(fit_seed)
    cv_perf <- .run_replicated_cv(
      args$std_data,
      cv_k = cv_k,
      cv_repl = cv_repl,
      metric = cv_metric,
      cv_est_fun = cv_est_fun,
      par_cluster = cl,
      handler_args = handler_args)

    data.frame(lambda = lambda, alpha = alpha,
               cvavg = rowMeans(cv_perf),
               cvse = if (cv_repl > 1L) { apply(cv_perf, 1, sd) } else { 0 })
  }

  cv_curves <- mapply(
    args$alpha, args$lambda, args$enpy_lambda_inds, other_starts,
    SIMPLIFY = FALSE, USE.NAMES = FALSE,
    FUN = cv_est_dispatch)

  cv_curves <- do.call(rbind, cv_curves)

  # if fit_all is not TRUE, compute only the fit for the "best" lambda
  if (!isTRUE(fit_ses)) {
    fit_lambda <- lapply(args$alpha, function (alpha) {
      rows <- which((cv_curves$alpha - alpha)^2 < .Machine$double.eps)

      lambda_inds <- vapply(
        fit_ses, FUN.VALUE = numeric(1L),
        FUN = function (se_fact) {
          which(.cv_se_selection(cv_curves$cvavg[rows],
                                 cv_curves$cvse[rows], se_fact) == 'se_fact')
        })
      unique(cv_curves$lambda[rows[lambda_inds]])
    })

    fit_enpy_lambda_inds <- lapply(fit_lambda, seq_along)

    fits <- .pense_internal_multi(args, lambda_list = fit_lambda,
                                  enpy_lambda_inds_list = fit_enpy_lambda_inds)
  }

  structure(list(
    call = match.call(expand.dots = TRUE),
    bdp = stable_bdp,
    lambda = args$lambda,
    alpha = vapply(fits, FUN.VALUE = numeric(1L),
                   FUN = `[[`, 'alpha', USE.NAMES = FALSE),
    cvres = cv_curves,
    cv_measure = cv_measure_str,
    cv_repl = cv_repl,
    metrics = lapply(fits, function (f) { attr(f$estimates, 'metrics') }),
    estimates = unlist(lapply(fits, `[[`, 'estimates'),
                       recursive = FALSE, use.names = FALSE)),
    class = c('pense', 'pense_cvfit'))
}

#' @description `adapense_cv()` is a convenience wrapper to compute adaptive
#'    PENSE estimates.
#'
#' @details
#' `adapense_cv()` is a convenience wrapper which performs 3 steps:
#'
#' 1. compute preliminary estimates via
#'    `pense_cv(..., alpha = alpha_preliminary)`,
#' 2. computes the penalty loadings from the estimate `beta` with best
#'    prediction performance by
#'    `adapense_loadings = 1 / abs(beta)^exponent`, and
#' 3. compute the adaptive PENSE estimates via
#'    `pense_cv(..., penalty_loadings = adapense_loadings)`.
#'
#' @param alpha_preliminary `alpha` parameter for the preliminary estimate.
#' @param exponent the exponent for computing the penalty loadings based on
#'     the preliminary estimate.
#'
#' @return a list-like object as returned by [pense_cv()] plus the following
#'    \describe{
#'      \item{`preliminary`}{the CV results for the preliminary estimate.}
#'      \item{`exponent`}{exponent used to compute the penalty loadings.}
#'      \item{`penalty_loadings`}{penalty loadings used for the
#'                                adaptive PENSE estimate.}
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

  adapense <- pense_cv(x, y, alpha = alpha, penalty_loadings = pen_loadings,
                       ...)
  adapense$call <- call
  adapense$exponent <- exponent
  adapense$preliminary <- prelim
  adapense$penalty_loadings <- pen_loadings
  class(adapense) <- c('adapense', class(adapense))
  return(adapense)
}

## Perform some final input adjustments and call the internal C++ code.
.pense_internal <- function(x, y, alpha, lambda, enpy_lambda_inds,
                            penalty_loadings = NULL,
                            pense_opts, enpy_opts, optional_args) {
  # Create penalties-list, without sorting the lambda sequence
  penalties <- lapply(lambda, function (l) { list(lambda = l, alpha = alpha) })

  if (!is.null(penalty_loadings)) {
    optional_args$pen_loadings <- penalty_loadings
  }

  .Call(C_pense_regression, x, y, penalties, enpy_lambda_inds, pense_opts,
        enpy_opts, optional_args)
}

## A wrapper to call `.pense_internal` for a list of `alpha` values,
## a list of corresponding `lambda` sequences and ENPY lambda indices.
## This wrapper filters the individual starts to only include the starts
## for the appropriate `alpha` value, and adds extra information to the
## returned fit.
.pense_internal_multi <- function (args, alpha_seq = args$alpha,
                                   lambda_list = args$lambda,
                                   enpy_lambda_inds_list = args$enpy_lambda_inds) {
  mapply(alpha_seq, lambda_list, enpy_lambda_inds_list,
         SIMPLIFY = FALSE, USE.NAMES = FALSE,
         FUN = function (alpha, lambda, enpy_lambda_inds) {
           # If there are other individual starts, only use the ones with
           # correct `alpha`
           if (length(args$optional_args$individual_starts) > 0L) {
             args$optional_args$individual_starts <- lapply(
               args$optional_args$individual_starts,
               FUN = .filter_list, what = 'alpha', value = alpha)
           }

           fit <- .pense_internal(x = args$std_data$x, y = args$std_data$y,
                                  alpha = alpha,
                                  lambda = lambda,
                                  enpy_lambda_inds = enpy_lambda_inds,
                                  penalty_loadings = args$penalty_loadings,
                                  pense_opts = args$pense_opts,
                                  enpy_opts = args$enpy_opts,
                                  optional_args = args$optional_args)

           # Flatten the list of estimates and un-standardize
           fit$estimates <- lapply(
             unlist(fit$estimates, recursive = FALSE),
             function (ests) {
               args$restore_coef_length(
                 args$std_data$unstandardize_coefs(ests))
             })

           # Handle metrics
           fit$estimates <- .metrics_attrib(fit$estimates, fit$metrics)
           fit$lambda <- unlist(vapply(fit$estimates, FUN.VALUE = numeric(1),
                                       FUN = `[[`, 'lambda'),
                                use.names = FALSE, recursive = FALSE)
           fit$alpha <- alpha
           fit
         })
}
