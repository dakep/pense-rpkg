#' Deprecated
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge('deprecated')}
#'
#' Options for computing initial estimates via ENPY.
#' Superseded by [enpy_options()].
#'
#' @details # Warning
#' Do not use this function in new code.
#' It may be removed from future versions of the package.
#'
#' @param keep_solutions how many initial estimates should be kept to perform
#'      full PENSE iterations?
#' @param psc_method The method to use for computing the principal sensitivity
#'      components. See details for the possible choices.
#' @param maxit maximum number of refinement iterations.
#' @param maxit_pense_refinement **ignored.** Maximum number of PENSE iterations to refine
#'      initial estimator.
#' @param eps **ignored.** Numeric tolerance for convergence.
#' @param psc_keep proportion of observations to keep based on the PSC scores.
#' @param resid_keep_method How to clean the data based on large residuals.
#'      If \code{"proportion"}, observations with the smallest
#'      \code{resid_keep_prop} residuals will be retained.
#'      If \code{"threshold"}, all observations with scaled residuals smaller
#'      than the threshold \code{resid_keep_thresh} will be retained.
#' @param resid_keep_prop,resid_keep_thresh proportion or threshold for
#'      observations to keep based on their residual.
#' @param mscale_eps,mscale_maxit **ignored.** Maximum number of iterations and numeric
#'      tolerance for the M-scale.
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

#' Deprecated
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge('deprecated')}
#'
#' Options for computing EN estimates.
#'
#' @details # Warning
#' Do not use these functions in new code.
#' They may be removed from future versions of the package.
#'
#' @name deprecated_en_options
#' @family deprecated functions
NULL


#' @param use_gram **ignored.** Should the Gram matrix be pre-computed.
#' @param eps **ignored.** Numeric tolerance for convergence.
#'
#' @export
#' @describeIn deprecated_en_options Superseded by [en_lars_options()].
#' @importFrom lifecycle deprecate_warn
en_options_aug_lars <- function (use_gram = c("auto", "yes", "no"), eps = 1e-12) {
  deprecate_warn('2.0.0', 'en_options_aug_lars()', with = 'en_lars_options()')
  en_lars_options()
}

#' @param maxit maximum number of iterations allowed.
#' @param eta_mult multiplier to increase eta at each iteration.
#' @param eta_start **ignored.** The start value for eta.
#' @param eta_start_numerator if \code{eta_start} is missing, it is defined
#'      by \code{eta_start = eta_start_numerator / lambda}.
#' @param preconditioner **ignored.** Preconditioner for the numerical solver. If none,
#'      a standard solver will be used, otherwise the faster preconditioned
#'      conjugate gradient is used.
#' @param verbosity **ignored.**
#'
#' @export
#' @describeIn deprecated_en_options Superseded by [en_dal_options()]
en_options_dal <- function (maxit = 100, eps = 1e-8, eta_mult = 2, eta_start_numerator = 1e-2,
                            eta_start, preconditioner = c("approx", "none", "diagonal"), verbosity = 0) {
  deprecate_warn('2.0.0', 'en_options_dal()', with = 'en_dal_options()')
  en_dal_options(max_it = .as(maxit[[1L]], 'integer'), eta_multiplier = .as(eta_mult[[1L]], 'numeric'),
                 eta_start_conservative = .as(eta_start_numerator[[1L]], 'numeric'),
                 eta_start_aggressive = .as(eta_start_numerator[[1L]], 'numeric'))
}


#' Deprecated
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge('deprecated')}
#'
#' Additional options for computing penalized EN S-estimates.
#' Superseded by [mm_algorithm_options()] and options supplied directly to [pense()].
#'
#' @details # Warning
#' Do not use this function in new code.
#' It may be removed from future versions of the package.
#'
#' @param delta desired breakdown point of the resulting estimator.
#' @param maxit maximum number of iterations allowed.
#' @param eps numeric tolerance for convergence.
#' @param mscale_eps,mscale_maxit maximum number of iterations and numeric
#'      tolerance for the M-scale.
#' @param cc **ignored.** Tuning constant for the S-estimator. Default is to chosen based
#'      on the breakdown point \code{delta}. Should never have to be changed.
#' @param en_correction **ignored.** Should the corrected EN estimator be used to choose
#'      the optimal lambda with CV.
#'      If \code{TRUE}, as by default, the estimator is "bias corrected".
#' @param verbosity **ignored.** Verbosity of the algorithm.
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

#' Deprecated
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge('deprecated')}
#'
#' Additional options for computing penalized EN MM-estimates.
#' Superseded by [mm_algorithm_options()] and options supplied directly to [pensem_cv()].
#'
#' @details # Warning
#' Do not use this function in new code.
#' It may be removed from future versions of the package.
#'
#' @param cc **ignored.** Tuning constant for the M-estimator.
#' @param maxit maximum number of iterations allowed.
#' @param eps **ignored.** Numeric tolerance for convergence.
#' @param adjust_bdp **ignored.** Should the breakdown point be adjusted based on the
#'      effective degrees of freedom?
#' @param en_correction **ignored.** Should the corrected EN estimator be used to choose
#'      the optimal lambda with CV.
#'      If \code{TRUE}, as by default, the estimator is "bias corrected".
#' @param verbosity **ignored.** Verbosity of the algorithm.
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

#' Deprecated
#'
#' @description
#' \Sexpr[results=rd, stage=render]{lifecycle::badge('deprecated')}
#'
#' Compute initial estimates for EN S-estimates using ENPY.
#' Superseded by [enpy_initial_estimates()].
#'
#' @details # Warning
#' Do not use this function in new code.
#' It may be removed from future versions of the package.
#'
#' @param x data matrix with predictors.
#' @param y response vector.
#' @param alpha,lambda EN penalty parameters (NOT adjusted for the number of
#'      observations in \code{x}).
#' @param delta desired breakdown point of the resulting estimator.
#' @param cc tuning constant for the S-estimator. Default is to chosen based
#'      on the breakdown point \code{delta}. Should never have to be changed.
#' @param options **ignored.** Additional options for the initial estimator.
#' @param en_options **ignored.** Additional options for the EN algorithm.
#'
#' @return \item{coeff}{A numeric matrix with one initial coefficient per column}
#'         \item{objF}{A vector of values of the objective function for the respective coefficient}
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
