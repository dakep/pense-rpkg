#' Deprecated Options for Initial Estimates
#'
#' \Sexpr[results=rd, stage=render]{lifecycle::badge('deprecated')}
#'
#' Superseded by [enpy_options()].
#'
#' @family deprecated functions
#'
#' @export
#'
#' @importFrom lifecycle deprecate_warn
initest_options <- function (keep_solutions = 5, psc_method = c("exact", "rr"), maxit = 10, maxit_pense_refinement = 5,
                             eps = 1e-6, psc_keep = 0.5, resid_keep_method = c("proportion", "threshold"),
                             resid_keep_prop = 0.6, resid_keep_thresh = 2, mscale_eps = 1e-8, mscale_maxit = 200) {
  deprecate_warn('2.0.0', 'initest_options()', with = 'enpy_options()')

  resid_keep_method <- match.arg(resid_keep_method)

  list(keepSolutions = .as(keep_solutions[[1L]], 'integer'), pscMethod = psc_method,
       maxit = .as(maxit[[1L]], 'integer'), keepResidualsMethod = resid_keep_method,
       keepResidualsThreshold = .as(resid_keep_thresh[[1L]], 'numeric'),
       keepResidualsProportion = .as(resid_keep_prop[[1L]], 'numeric'),
       keepPSCProportion = .as(psc_keep[[1L]], 'numeric'))
}

#' Deprecated Options for the EN Algorithms
#'
#' \Sexpr[results=rd, stage=render]{lifecycle::badge('deprecated')}
#'
#' @name deprecated_en_options
#' @family deprecated functions
NULL


#' @export
#' @describeIn deprecated_en_options Superseded by [en_lars_options()].
#' @importFrom lifecycle deprecate_warn
en_options_aug_lars <- function (use_gram = c("auto", "yes", "no"), eps = 1e-12) {
  deprecate_warn('2.0.0', 'en_options_aug_lars()', with = 'en_lars_options()')
  en_lars_options()
}

#' @export
#' @describeIn deprecated_en_options Superseded by [en_dal_options()]
en_options_dal <- function (maxit = 100, eps = 1e-8, eta_mult = 2, eta_start_numerator = 1e-2,
                            eta_start, preconditioner = c("approx", "none", "diagonal"), verbosity = 0) {
  deprecate_warn('2.0.0', 'en_options_dal()', with = 'en_dal_options()')
  en_dal_options(max_it = .as(maxit[[1L]], 'integer'), eta_multiplier = .as(eta_mult[[1L]], 'numeric'),
                 eta_start_conservative = .as(eta_start_numerator[[1L]], 'numeric'),
                 eta_start_aggressive = .as(eta_start_numerator[[1L]], 'numeric'))
}


#' Deprecated Additional Options for the Penalized EN S-estimator
#'
#' \Sexpr[results=rd, stage=render]{lifecycle::badge('deprecated')}
#'
#' Superseded by [mm_algorithm_options()] and options supplied directly to [pense()].
#'
#' @family deprecated functions
#'
#' @export
#'
#' @importFrom lifecycle deprecate_warn
pense_options <- function (delta = 0.25, maxit = 1000, eps = 1e-6, mscale_eps = 1e-8, mscale_maxit = 200,
                           verbosity = 0, cc = NULL, en_correction = TRUE) {
  deprecate_warn('2.0.0', 'pense_options()')
  list(maxit = .as(maxit[[1L]], 'integer'), bdp = .as(delta[[1L]], 'numeric'), eps = .as(eps[[1L]], 'numeric'),
       mscale_opts = mscale_algorithm_options(max_it = maxit, eps = mscale_eps))
}

#' Deprecated Additional Options for the Penalized EN MM-estimator
#'
#' \Sexpr[results=rd, stage=render]{lifecycle::badge('deprecated')}
#'
#' Superseded by [mm_algorithm_options()] and options supplied directly to [pense()].
#'
#' @family deprecated functions
#'
#' @export
#'
#' @importFrom lifecycle deprecate_warn
mstep_options <- function (cc = 3.44, maxit = 1000, eps = 1e-6, adjust_bdp = FALSE, verbosity = 0,
                           en_correction = TRUE) {
  deprecate_warn('2.0.0', 'mstep_options()')
  list(maxit = .as(maxit[[1L]], 'integer'))
}

#' Deprecated: ENPY Initial Estimates for EN S-Estimators
#'
#' \Sexpr[results=rd, stage=render]{lifecycle::badge('deprecated')}
#'
#' Superseded by [enpy_initial_estimates()].
#'
#' @family deprecated functions
#'
#' @export
#' @importFrom lifecycle deprecate_warn
enpy <- function(x, y, alpha, lambda, delta, cc, options, en_options) {
  deprecate_warn('2.0.0', 'enpy()', with = 'enpy_initial_estimates()')

  new_enpy <- enpy_initial_estimates(x, y, alpha = alpha, lambda = lambda, bdp = delta, cc = cc)
  coefs <- unlist(lapply(new_enpy, function (est) {
    c(est$intercept, as.numeric(est$beta))
  }), recursive = FALSE, use.names = FALSE)
  coefs <- matrix(coefs, ncol = length(new_enpy))
  objfv <- unlist(lapply(new_enpy, `[[`, 'objf_value'), recursive = FALSE, use.names = FALSE)
  list(coeff = coefs, objf = objfv)
}
