
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
                          en_algorithm_opts, keep_residuals_measure = c('threshold', 'proportion'),
                          keep_residuals_proportion = 0.5,  keep_residuals_threshold = 2,
                          retain_best_factor = 2, retain_max = 500) {
  list(max_it = .as(max_it[[1L]], 'integer'),
       en_options = if (missing(en_algorithm_opts)) { NULL } else { en_algorithm_opts },
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
.full_mscale_algo_options <- function (bdp, cc = NULL, mscale_opts) {
  mscale_opts$delta <- .as(bdp, 'numeric')
  if (isTRUE(mscale_opts$delta < .Machine$double.eps) || isTRUE(mscale_opts$delta > 0.5)) {
    stop("`bdp` is outside of 0 and 0.5")
  }

  mscale_opts$cc <- if (is.null(cc)) { .bisquare_consistency_const(mscale_opts$delta) } else { .as(cc, 'numeric') }
  return(mscale_opts)
}

#' MM-Algorithm to Compute Penalized Elastic Net S- and M-Estimates
#'
#' Additional options for the MM algorithm to compute EN S- and M-estimates.
#'
#' @param max_it maximum number of iterations.
#' @param tightening how to make inner iterations more precise as the algorithm approaches a local minimum.
#' @param tightening_steps for *adaptive* tightening strategy, how often to tighten until the desired tolerance
#'    is attained.
#' @param en_algorithm_opts options for the inner LS-EN algorithm. See [en_algorithm_options] for details.
#'
#' @return options for the MM algorithm.
#' @export
mm_algorithm_options <- function (max_it = 500, tightening = c('adaptive', 'exponential', 'none'),
                                  tightening_steps = 10, en_algorithm_opts) {
  list(algorithm = 'mm', max_it = .as(max_it[[1L]], 'integer'),
       tightening = .tightening_id(match.arg(tightening)),
       tightening_steps = .as(tightening_steps[[1L]], 'integer'),
       en_options = if (missing(en_algorithm_opts)) { NULL } else { en_algorithm_opts })
}


#' Control the Algorithm to Compute (Weighted) Least-Squares Elastic Net Estimates
#'
#' The package supports different algorithms to compute the EN estimate for weighted LS loss functions.
#' Each algorithm has certain characteristics that make it useful for some problems.
#' To select a specific algorithm and adjust the options, use any of the `en_***_options` functions.
#'
#' * [en_lars_options()]: Use the tuning-free LARS algorithm. This computes _exact_ (up to numerical errors) solutions
#'    to the EN-LS problem. It is not iterative and therefore can not benefit from approximate
#'    solutions, but in turn guarantees that a solution will be found.
#' * [en_admm_options()]: Use an iterative ADMM-type algorithm which needs \eqn{O(n p)} operations per iteration and
#'    converges sub-linearly.
#' * [en_dal_options()]: Use the iterative Dual Augmented Lagrangian (DAL) method. DAL needs \eqn{O(n^3 p^2)} operations
#'    per iteration, but converges exponentially.
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
  tau <- if (missing(step_size) || is.null(step_size)) { -1 } else { .as(step_size[[1L]], 'numeric') }
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
       eta_start_numerator_conservative = .as(eta_start_conservative[[1L]], 'numeric'),
       eta_start_numerator_aggressive = .as(eta_start_aggressive[[1L]], 'numeric'),
       lambda_relchange_aggressive = .as(lambda_relchange_aggressive[[1L]], 'numeric'),
       eta_multiplier = .as(eta_multiplier[[1L]], 'numeric'))
}

#' Ridge optimizer using an Augmented data matrix.
#' Only available for Ridge problems (`alpha=0``) and selected automatically in this case.
#'
#' @keywords internal
en_ridge_options <- function () {
  list(algorithm = 'augridge', sparse = FALSE)
}

## Check if the selected EN algorithm can handle the given penalty.
.check_en_algorithm <- function (en_algorithm_opts, alpha) {
  if (en_algorithm_opts$algorithm == 'dal') {
    if(!isTRUE(alpha > 0)) {
      warning('The DAL algorithm can not handle a Ridge penalty. Using default
              algorithm as fallback.')
      return(FALSE)
    }
  }
  return(TRUE)
}

## Choose the appropriate ADMM algorithm type based on the penalty and
## the problem size.
.choose_admm_algorithm <- function (en_algorithm_opts, alpha, x) {
  return(en_algorithm_opts)
}

.en_algorithm_id <- function (en_algorithm_opts) {
  switch (en_algorithm_opts$algorithm,
          admm = switch(en_algorithm_opts$admm_type, `var-stepsize` = 2L, 1L),
          dal = 3L,
          augridge = 4L,
          lars = 5L,
          5L)
}

.pense_algorithm_id <- function (opts) {
  switch (opts$algorithm,
          admm = 2L,
          mm = 1L,
          1L)
}

.regmest_algorithm_id <- function (opts) {
  switch (opts$algorithm,
          mm = 1L,
          1L)
}

.tightening_id <- function (tightening) {
  switch (tightening, exponential = 1L, adaptive = 2L, 0L)
}

#' @importFrom rlang warn
## also adds an element `sparse` to the returned list, which is the required sparsity parameter!
.select_en_algorithm <- function (en_options, alpha, sparse, eps) {
  if (!is.null(en_options)) {
    # Check if the EN algorithm can handle the given `alpha`
    if (en_options$algorithm == 'dal' && !isTRUE(alpha > 0)) {
      warn("The DAL algorithm can not handle a Ridge penalty. Using default algorithm as fallback.")
      en_options <- NULL
    } else if (en_options$algorithm == 'ridge' && !isTRUE(alpha < .Machine$double.eps)) {
      warn("The Ridge algorithm can only handle a Ridge penalty. Using default algorithm as fallback.")
      en_options <- NULL
    }
  }

  if (is.null(en_options)) {
    en_options <- if (alpha > 0) {
      en_lars_options()
    } else {
      en_ridge_options()
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
