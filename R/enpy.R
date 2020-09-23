#' ENPY Initial Estimates for EN S-Estimators
#'
#' Compute initial estimates for the EN S-estimator using the EN-PY procedure.
#'
#' If these manually computed initial estimates are intended as starting points for [pense()], they are by default
#' *shared* for all penalization levels.
#' To restrict the use of the initial estimates to the penalty level they were computed for, use
#' `as_starting_point(..., specific = TRUE)`. See [as_starting_point()] for details.
#'
#' @param y vector of response values of length `n`.
#' @inheritParams pense
#' @param enpy_opts options for the EN-PY algorithm, created with the [enpy_options()] function.
#' @param lambda a vector of positive values of penalization levels.
#' @param cc cutoff value for the bisquare rho function. By default, chosen to yield a consistent estimate for the
#'    Normal distribution.
#'
#' @family functions for initial estimates
#'
#' @export
#'
#' @references Cohen Freue, G.V.; Kepplinger, D.; Salibián-Barrera, M.; Smucler, E.
#'    Robust elastic net estimators for variable selection and identification of proteomic biomarkers.
#'    *Ann. Appl. Stat.* **13** (2019), no. 4, 2065--2090 \doi{10.1214/19-AOAS1269}
#'
#'
#' @importFrom lifecycle deprecate_warn deprecate_stop deprecated is_present
#' @importFrom rlang abort
enpy_initial_estimates <- function (x, y, alpha, lambda, bdp = 0.25, cc, intercept = TRUE, penalty_loadings,
                                    enpy_opts = enpy_options(), mscale_opts = mscale_algorithm_options(),
                                    eps = 1e-6, sparse = FALSE, ncores = 1L) {
  y <- .as(y, 'numeric')
  x_dim <- dim(x)

  if (length(y) != x_dim[[1L]]) {
    abort("Number of observations in `x` and `y` does not match.")
  } else if (x_dim[[2L]] <= 1L) {
    abort("`x` must be a matrix with at least 1 column.")
  }

  if (missing(cc)) {
    cc <- NULL
  }

  alpha <- .as(alpha[[1L]], 'numeric')
  if (alpha < 0 || alpha > 1) {
    abort("`alpha` is outside 0 and 1.")
  } else if (alpha < sqrt(.Machine$double.eps)) {
    alpha <- 0
  }

  sparse <- .as(sparse[[1L]], 'logical')

  # Check EN algorithm for ENPY
  enpy_opts$num_threads <- max(1L, .as(ncores[[1L]], 'integer'))
  enpy_opts$eps <- .as(eps[[1L]], 'numeric')
  enpy_opts$en_options <- .select_en_algorithm(enpy_opts$en_options, alpha, sparse, enpy_opts$eps)
  sparse <- enpy_opts$en_options$sparse

  if (enpy_opts$num_threads > 1L && !isTRUE(.k_multithreading_support)) {
    warn("Multithreading not supported. Using only 1 core.")
  }

  optional_args <- list()
  if (!missing(penalty_loadings) && !is.null(penalty_loadings)) {
    checked_pls <- .prepare_penalty_loadings(penalty_loadings, x, alpha, sparse = sparse, stop_all_infinite = TRUE)
    optional_args$penalty_loadings <- checked_pls$loadings
    restore_coef_length <- checked_pls$restore_fun
    x <- checked_pls$trimmed_x
  } else {
    restore_coef_length <- function (coef) coef
    optional_args$penalty_loadings <- NULL
  }

  lambda <- sort(.as(lambda, 'numeric'), decreasing = TRUE)
  penalties <- lapply(lambda, function (lambda) { list(alpha = alpha, lambda = lambda) })
  s_loss_params <- list(mscale = .full_mscale_algo_options(bdp = bdp, cc = cc, mscale_opts = mscale_opts),
                        intercept = intercept)

  res <- .Call(C_penpy, x, drop(y), penalties, s_loss_params, enpy_opts, optional_args)

  res <- lapply(res, function (lambda_res) {
    lapply(.metrics_attrib(lambda_res$estimates, lambda_res$metrics), function (est) {
      structure(restore_coef_length(est), class = c('shared_starting_point', 'starting_point'))
    })
  })

  structure(unlist(res, recursive = FALSE, use.names = FALSE), class = c('enpy_starting_points', 'starting_points'))
}

#' Principal Sensitivity Components
#'
#' Compute Principal Sensitivity Components for Elastic Net Regression
#'
#' @param y vector of response values of length `n`.
#' @inheritParams pense
#' @param en_algorithm_opts options for the LS-EN algorithm. See [en_algorithm_options] for details.
#' @param method defunct. PSCs are always computed for EN estimates. For the PY procedure for unpenalized estimation
#'    use package [pyinit](https://cran.r-project.org/package=pyinit).
#'
#' @return a list of principal sensitivity components, one per element in `lambda`. Each PSC is itself a list
#'    with items `lambda`, `alpha`, and `pscs`.
#'
#' @family functions for initial estimates
#'
#' @export
#'
#' @references Cohen Freue, G.V.; Kepplinger, D.; Salibián-Barrera, M.; Smucler, E.
#'    Robust elastic net estimators for variable selection and identification of proteomic biomarkers.
#'    *Ann. Appl. Stat.* **13** (2019), no. 4, 2065--2090 \doi{10.1214/19-AOAS1269}
#' @references Pena, D., and Yohai, V.J.
#'    A Fast Procedure for Outlier Diagnostics in Large Regression Problems.
#'    *J. Amer. Statist. Assoc.* **94** (1999). no. 446, 434--445. \doi{10.2307/2670164}
#'
#' @importFrom lifecycle deprecate_warn deprecate_stop deprecated is_present
#' @importFrom rlang abort
prinsens <- function (x, y, alpha, lambda, intercept = TRUE, penalty_loadings, en_algorithm_opts,
                      eps = 1e-6, sparse = FALSE, ncores = 1L, method = deprecated()) {
  y <- .as(y, 'numeric')
  x_dim <- dim(x)

  if (length(y) != x_dim[[1L]]) {
    abort("Number of observations in `x` and `y` does not match.")
  } else if (x_dim[[2L]] <= 1L) {
    abort("`x` must be a matrix with at least 1 column.")
  }

  if (is_present(method)) {
    if (method == 'en') {
      deprecate_warn('2.0.0', 'prinsens(method = )')
    } else {
      deprecate_stop('2.0.0', 'prinsens(method = )', details = 'For unpenalized estimates use the pyinit package.')
    }
  }

  sparse <- .as(sparse[[1L]], 'logical')

  if (missing(en_algorithm_opts)) {
    en_algorithm_opts <- NULL
  }

  alpha <- .as(alpha[[1L]], 'numeric')
  if (alpha < 0 || alpha > 1) {
    abort("`alpha` is outside 0 and 1.")
  } else if (alpha < sqrt(.Machine$double.eps)) {
    alpha <- 0
  }

  en_algorithm_opts <- .select_en_algorithm(en_algorithm_opts, alpha, sparse, eps)
  sparse <- en_algorithm_opts$sparse

  optional_args <- list(intercept = .as(intercept[[1L]], 'logical'),
                        num_threads = max(1L, .as(ncores[[1L]], 'integer')))

  if (!missing(penalty_loadings) && !is.null(penalty_loadings)) {
    checked_pls <- .prepare_penalty_loadings(penalty_loadings, x, alpha, sparse = sparse, stop_all_infinite = TRUE)
    optional_args$penalty_loadings <- checked_pls$loadings
    x <- checked_pls$trimmed_x
  } else {
    optional_args$penalty_loadings <- NULL
  }

  lambda <- sort(.as(lambda, 'numeric'), decreasing = TRUE)
  penalties <- lapply(lambda, function (lambda) { list(alpha = alpha, lambda = lambda) })
  eps <- .as(eps[[1L]], 'numeric')

  if (optional_args$num_threads > 1L && !isTRUE(.k_multithreading_support)) {
    warn("Multithreading not supported. Using only 1 core.")
  }

  res <- .Call(C_pscs, x, drop(y), penalties, en_algorithm_opts, optional_args)

  mapply(penalties, res, FUN = function (pen, pscs) {
    pen$pscs <- pscs
    return(pen)
  }, SIMPLIFY = FALSE, USE.NAMES = FALSE)
}
