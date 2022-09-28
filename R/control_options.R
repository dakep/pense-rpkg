#' Options for the ENPY Algorithm
#'
#' Additional control options for the elastic net Peña-Yohai procedure.
#'
#' The EN-PY procedure for computing initial estimates iteratively cleans the data of observations with possibly
#' outlying residual or high leverage. Least-squares elastic net (LS-EN) estimates are computed on the possibly clean
#' subsets. At each iteration, the Principal Sensitivity Components are computed to remove observations with potentially
#' high leverage. Among all the LS-EN estimates, the estimate with smallest M-scale of the residuals is selected.
#' Observations with largest residual for the selected estimate are removed and the next iteration is started.
#'
#' @param max_it maximum number of EN-PY iterations.
#' @param en_algorithm_opts options for the LS-EN algorithm. See [en_algorithm_options] for details.
#' @param keep_psc_proportion how many observations should to keep based on the Principal Sensitivity Components.
#' @param keep_residuals_measure how to determine what observations to keep, based on their residuals.
#'    If `proportion`, a fixed number of observations is kept.
#'    If `threshold`, only observations with residuals below the threshold are kept.
#' @param keep_residuals_proportion proportion of observations to kept based on their residuals.
#' @param keep_residuals_threshold only observations with (standardized) residuals less than this threshold are kept.
#' @param retain_best_factor only keep candidates that are within this factor of the best candidate. If `<= 1`, only
#'    keep candidates from the last iteration.
#' @param retain_max maximum number of candidates, i.e., only the best `retain_max` candidates are retained.
#'
#' @return options for the ENPY algorithm.
#' @export
enpy_options <- function (max_it = 10, keep_psc_proportion = 0.5,
                          en_algorithm_opts,
                          keep_residuals_measure = c('threshold', 'proportion'),
                          keep_residuals_proportion = 0.5,
                          keep_residuals_threshold = 2,
                          retain_best_factor = 2, retain_max = 500) {
  list(max_it = .as(max_it[[1L]], 'integer'),
       en_options = if (missing(en_algorithm_opts)) {
         NULL
       } else {
         en_algorithm_opts
       },
       keep_psc_proportion = .as(keep_psc_proportion[[1L]], 'numeric'),
       use_residual_threshold = match.arg(keep_residuals_measure) == 'threshold',
       keep_residuals_proportion = .as(keep_residuals_proportion[[1L]], 'numeric'),
       keep_residuals_threshold = .as(keep_residuals_threshold[[1L]], 'numeric'),
       retain_best_factor = .as(retain_best_factor[[1L]], 'numeric'),
       retain_max = .as(retain_max[[1L]], 'integer'))
}

#' Options for the M-scale Estimation Algorithm
#'
#' @param max_it maximum number of iterations.
#' @param eps numerical tolerance to check for convergence.
#'
#' @return options for the M-scale estimation algorithm.
#' @export
mscale_algorithm_options <- function (max_it = 200, eps = 1e-8) {
  list(max_it = .as(max_it[[1L]], 'integer'),
       eps = .as(eps[[1L]], 'numeric'))
}


## Full options for the M-scale Estimation Algorithm
##
## @param mscale_opts public control options created by [mscale_algorithm_options].
## @param bdp the breakdown point, i.e., `delta` in the M-estimation equation.
## @param cc the cutoff threshold for the bisquare rho function.
##
## @return full options for the M-scale estimation algorithm.
#' @importFrom rlang is_missing
.full_mscale_algo_options <- function (bdp, cc = NULL, mscale_opts) {
  mscale_opts$delta <- .as(bdp, 'numeric')
  if (isTRUE(mscale_opts$delta < 0.05) ||
      isTRUE(mscale_opts$delta > 0.5)) {
    stop("`bdp` is outside of 0.05 and 0.5")
  }

  mscale_opts$cc <- if (is_missing(cc) || is.null(cc)) {
    .bisquare_consistency_const(mscale_opts$delta)
  } else {
    .as(cc, 'numeric')
  }
  return(mscale_opts)
}

#' MM-Algorithm to Compute Penalized Elastic Net S- and M-Estimates
#'
#' Additional options for the MM algorithm to compute EN S- and M-estimates.
#'
#' @param max_it maximum number of iterations.
#' @param tightening how to make inner iterations more precise as the algorithm
#'   approaches a local minimum.
#' @param tightening_steps for *adaptive* tightening strategy, how often to
#'   tighten until the desired tolerance is attained.
#' @param en_algorithm_opts options for the inner LS-EN algorithm.
#'   See [en_algorithm_options] for details.
#'
#' @return options for the MM algorithm.
#' @seealso cd_algorithm_options for a direct optimization of the non-convex
#'   PENSE loss.
#' @export
mm_algorithm_options <- function (max_it = 500,
                                  tightening = c('adaptive', 'exponential',
                                                 'none'),
                                  tightening_steps = 2,
                                  en_algorithm_opts) {
  list(algorithm = 'mm', max_it = .as(max_it[[1L]], 'integer'),
       tightening = .tightening_id(match.arg(tightening)),
       tightening_steps = .as(tightening_steps[[1L]], 'integer'),
       en_options = if (missing(en_algorithm_opts)) {
         NULL
       } else {
         en_algorithm_opts
       })
}

#' Coordinate Descent (CD) Algorithm to Compute Penalized Elastic Net
#' S-estimates
#'
#' Set options for the CD algorithm to compute adaptive EN S-estimates.
#'
#' @param max_it maximum number of iterations.
#' @param reset_it number of iterations after which the residuals are
#'   re-computed from scratch, to prevent numerical drifts from incremental
#'   updates.
#' @param linesearch_steps maximum number of steps used for line search.
#' @param linesearch_mult multiplier to adjust the step size in the line
#'   search.
#'
#' @return options for the CD algorithm to compute (adaptive) PENSE estimates.
#' @seealso mm_algorithm_options to optimize the non-convex PENSE objective
#'   function via a sequence of convex problems.
#' @export
#' @importFrom rlang abort
cd_algorithm_options <- function (max_it = 1000, reset_it = 8,
                                  linesearch_steps = 4,
                                  linesearch_mult = 0.5) {
  opts <- list(algorithm = 'cd',
               max_it = .as(max_it[[1L]], 'integer'),
               linesearch_steps = .as(linesearch_steps[[1L]], 'integer'),
               linesearch_mult = .as(linesearch_mult[[1L]], 'numeric'),
               reset_it = .as(reset_it[[1L]], 'integer'))

  if (opts$linesearch_mult <= 0 || opts$linesearch_mult >= 1) {
    abort("`linesearch_mult` must be between 0 and 1.")
  }
  if (opts$linesearch_steps < 1L) {
    abort("`linesearch_steps` must be at least 1.")
  }
  if (opts$max_it < 2L) {
    abort("`max_it` must be greater than 1.")
  }
  if (opts$reset_it < 2L) {
    abort("`reset_it` must be greater than 1.")
  }

  opts
}

#' Control the Algorithm to Compute (Weighted) Least-Squares Elastic
#' Net Estimates
#'
#' The package supports different algorithms to compute the EN estimate
#' for weighted LS loss functions.
#' Each algorithm has certain characteristics that make it useful for some
#' problems.
#' To select a specific algorithm and adjust the options, use any of
#' the `en_***_options` functions.
#'
#' * [en_lars_options()]: Use the tuning-free LARS algorithm.
#'    This computes _exact_ (up to numerical errors) solutions to the EN-LS
#'    problem.
#'    It is not iterative and therefore can not benefit from approximate
#'    solutions, but in turn guarantees that a solution will be found.
#' * [en_cd_options()]: Use an iterative coordinate descent algorithm which
#'    needs \eqn{O(n p)} operations per iteration and converges sub-linearly.
#' * [en_admm_options()]: Use an iterative ADMM-type algorithm which needs
#'    \eqn{O(n p)} operations per iteration and converges sub-linearly.
#' * [en_dal_options()]: Use the iterative Dual Augmented Lagrangian (DAL)
#'    method.
#'    DAL needs \eqn{O(n^3 p^2)} operations per iteration, but converges
#'    exponentially.
#'
#' @name en_algorithm_options
NULL

#' Use the LARS Elastic Net Algorithm
#'
#' @family EN algorithms
#' @export
en_lars_options <- function () {
  list(algorithm = 'lars')
}

#' Use Coordinate Descent to Solve Elastic Net Problems
#'
#' @param max_it maximum number of iterations.
#' @param reset_it number of iterations after which the residuals are
#'   re-computed from scratch, to prevent numerical drifts from incremental
#'   updates.
#' @family EN algorithms
#' @export
en_cd_options <- function (max_it = 1000, reset_it = 8) {
  list(algorithm = 'cdls',
       max_it = .as(max_it[[1L]], 'integer'),
       reset_it = .as(reset_it[[1L]], 'integer'))
}

#' Use the ADMM Elastic Net Algorithm
#'
#'
#' @param max_it maximum number of iterations.
#' @param step_size step size for the algorithm.
#' @param acceleration acceleration factor for linearized ADMM.
#'
#' @return options for the ADMM EN algorithm.
#' @family EN algorithms
#' @export
en_admm_options <- function (max_it = 1000, step_size, acceleration = 1) {
  tau <- if (missing(step_size) || is.null(step_size)) {
    -1
  } else {
    .as(step_size[[1L]], 'numeric')
  }
  list(algorithm = 'admm',
       admm_type = 'linearized',
       max_it = .as(max_it[[1L]], 'integer'),
       accelerate = .as(acceleration[[1L]], 'numeric'),
       prox_opts = list(tau = tau), tau = tau)
}

#' Use the DAL Elastic Net Algorithm
#'
#' @param max_it maximum number of (outer) iterations.
#' @param max_inner_it maximum number of (inner) iterations in each outer iteration.
#' @param eta_multiplier multiplier for the barrier parameter. In each iteration, the barrier must be more restrictive
#'    (i.e., the multiplier must be > 1).
#' @param eta_start_conservative conservative initial barrier parameter. This is used if the previous penalty is
#'    undefined or too far away.
#' @param eta_start_aggressive aggressive initial barrier parameter. This is used if the previous penalty is close.
#' @param lambda_relchange_aggressive how close must the lambda parameter from the previous penalty term be to use
#'    an aggressive initial barrier parameter (i.e., what constitutes "too far").
#'
#' @return options for the DAL EN algorithm.
#' @family EN algorithms
#' @export
en_dal_options <- function (max_it = 100, max_inner_it = 100, eta_multiplier = 2,
                            eta_start_conservative = 0.01, eta_start_aggressive = 1,
                            lambda_relchange_aggressive = 0.25) {
  list(algorithm = 'dal',
       sparse = TRUE,
       max_it = .as(max_it[[1L]], 'integer'),
       max_inner_it = .as(max_inner_it[[1L]], 'integer'),
       eta_start_numerator_conservative = .as(eta_start_conservative[[1L]],
                                              'numeric'),
       eta_start_numerator_aggressive = .as(eta_start_aggressive[[1L]],
                                            'numeric'),
       lambda_relchange_aggressive = .as(lambda_relchange_aggressive[[1L]],
                                         'numeric'),
       eta_multiplier = .as(eta_multiplier[[1L]], 'numeric'))
}

#' Ridge optimizer using an Augmented data matrix.
#' Only available for Ridge problems (`alpha=0``) and selected automatically in this case.
#'
#' @keywords internal
en_ridge_options <- function () {
  list(algorithm = 'augridge', sparse = FALSE)
}

## Choose the appropriate ADMM algorithm type based on the penalty and
## the problem size.
.choose_admm_algorithm <- function (en_algorithm_opts, alpha, x) {
  return(en_algorithm_opts)
}

.en_algorithm_id <- function (en_algorithm_opts) {
  switch (en_algorithm_opts$algorithm,
          admm = switch(en_algorithm_opts$admm_type,
                        `var-stepsize` = .k_en_algo_admm_var,
                        .k_en_algo_admm_lin),
          dal = .k_en_algo_dal,
          augridge = .k_en_algo_augridge,
          lars = .k_en_algo_lars,
          cdls = .k_en_algo_cdls,
          .k_en_algo_lars)
}

.pense_algorithm_id <- function (opts) {
  switch (opts$algorithm,
          cd = .k_pense_algo_cd,
          admm = .k_pense_algo_admm,
          mm = .k_pense_algo_mm,
          .k_pense_algo_mm)
}

.regmest_algorithm_id <- function (opts) {
  switch (opts$algorithm,
          mm = .k_regm_algo_mm,
          .k_regm_algo_mm)
}

.tightening_id <- function (tightening) {
  switch (tightening, exponential = 1L, adaptive = 2L, 0L)
}

#' @importFrom rlang warn
## also adds an element `sparse` to the returned list, which is the required sparsity parameter!
.select_en_algorithm <- function (en_options, alpha, sparse, eps) {
  if (!is.null(en_options)) {
    # Check if the EN algorithm can handle the given `alpha`
    if (!all(alpha > 0) && (identical(en_options$algorithm, 'dal') ||
                            identical(en_options$algorithm, 'cdls'))) {
      warn(paste("The chosen algorithm can not handle an EN penalty with alpha=0.",
                 "Using default algorithm as fallback."))
      en_options <- NULL
    } else if (identical(en_options$algorithm, 'augridge') && !all(alpha < .Machine$double.eps)) {
      warn(paste("The Ridge algorithm can only handle a Ridge penalty.",
                 "Using default algorithm as fallback."))
      en_options <- NULL
    } else if (identical(en_options$algorithm, 'lars') && all(alpha < .Machine$double.eps)) {
      en_options <- en_ridge_options()
    }
  }

  if (is.null(en_options)) {
    en_options <- if (all(alpha > 0)) {
      en_lars_options()
    } else if (all(alpha < .Machine$double.eps)) {
      en_ridge_options()
    } else {
      abort("Cannot compute PENSE for alpha=0 and alpha > 0 simultaneously.")
    }
  }

  # Set the numerical tolerance
  en_options$eps <- .as(eps[[1L]], 'numeric')

  if (en_options$algorithm == 'admm') {
    # Numerical tolerance for convergence must be tightened for ADMM-type algorithms.
    en_options$eps <- 0.1 * en_options$eps
  }

  en_options$algorithm <- .en_algorithm_id(en_options)
  # Set sparsity if not pre-determined
  if (is.null(en_options$sparse)) {
    en_options$sparse <- isTRUE(sparse[[1L]])
  }
  return(en_options)
}
