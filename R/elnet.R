#' Compute the Least Squares (Adaptive) Elastic Net Regularization Path
#'
#' Compute least squares EN estimates for linear regression with optional observation
#' weights and penalty loadings.
#'
#' The elastic net estimator for the linear regression model solves
#' the optimization problem
#'
#' \deqn{argmin_{\mu, \beta}
#'   (1/2n) \sum_i w_i (y_i - \mu - x_i' \beta)^2 +
#'   \lambda \sum_j 0.5 (1 - \alpha) \beta_j^2 + \alpha l_j |\beta_j|  }
#'
#' with observation weights \eqn{w_i} and penalty loadings \eqn{l_j}.
#'
#' @param x `n` by `p` matrix of numeric predictors.
#' @param y vector of response values of length `n`.
#'          For binary classification, `y` should be a factor with 2 levels.
#' @param alpha elastic net penalty mixing parameter with \eqn{0 \le \alpha \le 1}. `alpha = 1` is the LASSO penalty,
#'    and `alpha = 0` the Ridge penalty.
#' @param nlambda number of penalization levels.
#' @param lambda_min_ratio Smallest value of the penalization level as a fraction of the largest level (i.e., the
#'    smallest value for which all coefficients are zero). The default depends on the sample
#'    size relative to the number of variables and `alpha`. If more observations than variables
#'    are available, the default is `1e-3 * alpha`, otherwise `1e-2 * alpha`.
#' @param lambda optional user-supplied sequence of penalization levels. If given and not `NULL`, `nlambda` and
#'    `lambda_min_ratio` are ignored.
#' @param penalty_loadings a vector of positive penalty loadings (a.k.a. weights) for different penalization of each
#'    coefficient.
#' @param standardize standardize variables to have unit variance. Coefficients are always returned in original scale.
#' @param weights a vector of positive observation weights.
#' @param intercept include an intercept in the model.
#' @param sparse use sparse coefficient vectors.
#' @param en_algorithm_opts options for the EN algorithm. See [en_algorithm_options]
#'    for details.
#' @param eps numerical tolerance.
#' @param xtest deprecated. Instead, extract coefficients with [coef.pense_fit()] and compute predictions manually.
#' @param options deprecated. Use `en_algorithm_opts` instead.
#' @param correction defunct. Correction for EN estimates is not supported anymore.
#'
#' @return a list-like object with the following items
#'    \describe{
#'      \item{`lambda`}{the sequence of penalization parameters.}
#'      \item{`estimates`}{a list of estimates. Each estimate contains the following information:
#'         \describe{
#'           \item{`intercept`}{intercept estimate.}
#'           \item{`beta`}{beta (slope) estimate.}
#'           \item{`lambda`}{penalization level at which the estimate is computed.}
#'           \item{`alpha`}{*alpha* hyper-parameter at which the estimate is computed.}
#'           \item{`statuscode`}{if `> 0` the algorithm experienced issues when computing the estimate.}
#'           \item{`status`}{optional status message from the algorithm.}
#'         }
#'      }
#'      \item{`call`}{the original call.}
#'      \item{`predictions`}{if `xtest` was given, a matrix of predicted values. Each column corresponds to the
#'                           predictions from the estimate at the `lambda` value at the same index.}
#'    }
#'
#' @example examples/ls_elnet.R
#'
#' @aliases adaelnet adaen
#' @family functions for computing non-robust estimates
#' @seealso [pense()] for an S-estimate of regression with elastic net penalty.
#' @seealso [coef.pense_fit()] for extracting coefficient estimates.
#' @seealso [plot.pense_fit()] for plotting the regularization path.
#'
#' @export
#' @importFrom lifecycle deprecated
#' @importFrom rlang exec
elnet <- function(x, y, alpha, nlambda = 100, lambda_min_ratio, lambda, penalty_loadings, weights, intercept = TRUE,
                  en_algorithm_opts, sparse = FALSE, eps = 1e-6, standardize = TRUE,
                  correction = deprecated(), xtest = deprecated(), options = deprecated()) {
  call <- match.call(expand.dots = FALSE)
  args <- as.list(call[-1L])
  args$standardize <- isTRUE(standardize)  # Ignore standardize = 'cv_only'!
  args <- do.call(.elnet_args, args, envir = parent.frame())

  res <- .elnet_internal(args$std_data$x, args$std_data$y, alpha = args$alpha, lambda = args$lambda,
                         penalty_loadings = args$penalty_loadings, weights = args$weights,
                         intercept = args$intercept, optional_args = args$optional_args)

  res$estimates <- lapply(res$estimates, function (est) {
    args$restore_coef_length(args$std_data$unstandardize_coefs(est))
  })

  predictions <- if (is_present(xtest)) {
    matrix(unlist(lapply(res$estimates, function (est) {
      as.numeric(xtest %*% est$beta) + est$intercept
    }), use.names = FALSE, recursive = FALSE), ncol = length(res$estimates))
  } else {
    NULL
  }

  structure(list(estimates = .metrics_attrib(res$estimates, res$metrics), call = call, predictions = predictions,
                 lambda = unlist(lapply(res$estimates, `[[`, 'lambda'), use.names = FALSE, recursive = FALSE)),
            class = c('en', 'pense_fit'))
}


#' Cross-validation for Least-Squares (Adaptive) Elastic Net Estimates
#'
#' Perform (repeated) K-fold cross-validation for [elnet()].
#'
#' @inheritParams elnet
#' @template cv_params
#' @inheritDotParams elnet
#' @param ncores deprecated and not used anymore.
#'
#' @seealso [elnet()] for computing the LS-EN regularization path without cross-validation.
#' @seealso [pense_cv()] for cross-validation of S-estimates of regression with elastic net penalty.
#'
#' @return a list with components:
#'    \describe{
#'      \item{`lambda`}{the sequence of penalization levels.}
#'      \item{`cvres`}{data frame of average cross-validated performance.}
#'      \item{`cv_replications`}{matrix of cross-validated performance metrics, one column per replication.
#'                               Rows are in the same order as in `cvres`.}
#'      \item{`call`}{the original call.}
#'      \item{`estimates`}{the estimates fitted on the full data. Same format as returned by [elnet()].}
#'    }
#'
#' @family functions for computing non-robust estimates
#' @seealso [coef.pense_cvfit()] for extracting coefficient estimates.
#' @seealso [plot.pense_cvfit()] for plotting the CV performance or the regularization path.
#' @example examples/ls_elnet.R
#'
#' @importFrom lifecycle deprecate_warn deprecated is_present
#' @importFrom stats sd
#' @export
elnet_cv <- function (x, y, lambda, cv_k, cv_repl = 1, cv_metric = c('rmspe', 'tau_size', 'mape', 'auroc'),
                      fit_all = TRUE, cl = NULL, ncores = deprecated(), ...) {
  call <- match.call(expand.dots = TRUE)
  args <- do.call(.elnet_args, as.list(call[-1L]), envir = parent.frame())
  cv_k <- .as(cv_k, 'integer')
  cv_repl <- .as(cv_repl, 'integer')

  if (cv_k < 2L) {
    abort("`cv_k` must be greater than 1.")
  }

  if (cv_repl < 1L) {
    abort("`cv_repl` must be greater than 0.")
  }

  if (is_present(ncores)) {
    deprecate_warn('2.0.0', 'elnet_cv(ncores=)', details = "Please specify the `parallel` cluster with argument `cl`.")
  }
  if (!is.null(cl) && !is(cl, 'cluster')) {
    abort("`cl` must be a valid `parallel` cluster.")
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

  cv_fun <- function (train_data, test_ind) {
    cv_fit <- .elnet_internal(train_data$x, train_data$y, alpha = args$alpha, lambda = args$lambda,
                              penalty_loadings = args$penalty_loadings, weights = args$weights[-test_ind],
                              intercept = args$intercept, optional_args = args$optional_args)
    cv_fit$estimates
  }

  cv_perf <- .run_replicated_cv(args$std_data, cv_k = cv_k, cv_repl = cv_repl, metric = cv_metric, cv_est_fun = cv_fun,
                                par_cluster = cl)

  cv_perf_df <- data.frame(lambda = args$lambda, cvavg = rowMeans(cv_perf), cvse = 0)

  if (cv_repl > 1L) {
    cv_perf_df$cvse <- apply(cv_perf, 1, sd)
  }

  fit_lambda <- if (isTRUE(fit_all)) {
    args$lambda
  } else {
    with(cv_perf_df, lambda[[which.min(cv_perf_df$cvavg)]])
  }

  fit <- .elnet_internal(args$std_data$x, args$std_data$y, alpha = args$alpha, lambda = fit_lambda,
                         penalty_loadings = args$penalty_loadings, weights = args$weights, intercept =args$intercept,
                         optional_args = args$optional_args)
  fit$estimates <- lapply(fit$estimates, function (est) {
    args$restore_coef_length(args$std_data$unstandardize_coefs(est))
  })
  return(structure(list(call = call, cvres = cv_perf_df, cv_replications = cv_perf, cv_measure = cv_measure_str,
                        lambda = args$lambda, estimates = .metrics_attrib(fit$estimates, fit$metrics)),
                   class = c('pense_en', 'pense_cvfit')))
}

## Perform some final input adjustments and call the internal C++ code.
.elnet_internal <- function(x, y, alpha, lambda, penalty_loadings, weights, intercept, optional_args) {
  # Create penalties-list, without sorting the lambda sequence
  penalties <- lapply(lambda, function (l) { list(lambda = l, alpha = alpha) })
  intercept <- isTRUE(intercept)

  if (!is.null(penalty_loadings)) {
    optional_args$pen_loadings <- penalty_loadings
  }

  if (!is.null(weights)) {
    optional_args$obs_weights <- weights
  }

  return(.Call(C_lsen_regression, x, y, penalties, intercept, optional_args))
}

## Check and parse user-supplied arguments
#' @importFrom lifecycle deprecate_warn deprecated is_present
#' @importFrom rlang warn abort
#' @importFrom stats runif
.elnet_args <- function (x, y, alpha, nlambda = 100, lambda_min_ratio, lambda, penalty_loadings, weights,
                         intercept = TRUE, en_algorithm_opts, sparse = FALSE, eps = 1e-6, standardize = TRUE,
                         correction = deprecated(), xtest = deprecated(), options = deprecated(), ...) {
  optional_args <- list()

  if (is_present(options)) {
    deprecate_warn('2.0.0', 'elnet(options=)', 'elnet(en_algorithm_opts=)')
  }

  if (is_present(correction)) {
    deprecate_warn('2.0.0', 'elnet(correction=)', 'coef()')
  }

  if (is_present(xtest)) {
    deprecate_warn('2.0.0', 'elnet(xtest=)')
  }

  # Normalize input
  response <- .validate_response(y)
  y <- response$values
  x_dim <- dim(x)

  if (length(y) != x_dim[[1L]]) {
    abort("Number of observations does not match between `x` and `y`.")
  } else if (x_dim[[2L]] <= 1L) {
    abort("`x` must be a matrix with at least 2 columns.")
  }

  intercept <- !isFALSE(intercept)
  standardize <- if (is.character(standardize)) {
    if (is.na(pmatch(standardize[[1L]], 'cv_only'))) {
      abort("`standardize` must be either TRUE/FALSE or \"cv_only\".")
    } else {
      standardize <- 'cv_only'
    }
  } else {
    isTRUE(standardize)
  }

  alpha <- .as(alpha[[1L]], 'numeric')
  if (alpha < 0 || alpha > 1) {
    abort("`alpha` is outside 0 and 1.")
  } else if (alpha < sqrt(.Machine$double.eps)) {
    alpha <- 0
  }

  sparse <- .as(sparse[[1L]], 'logical')

  # Check EN algorithm
  if (missing(en_algorithm_opts)) {
    en_algorithm_opts <- NULL
  }
  optional_args$en_options <- .select_en_algorithm(en_algorithm_opts, alpha, sparse, eps)
  sparse <- optional_args$en_options$sparse

  weights <- if (!missing(weights)) {
    if (length(weights) != x_dim[[1L]]) {
      abort("Observation weights are not the same length as `y`.")
    }
    .as(weights, 'numeric')
  } else {
    NULL
  }

  # Check penalty loadings
  if (!missing(penalty_loadings) && !is.null(penalty_loadings)) {
    checked_pls <- .prepare_penalty_loadings(penalty_loadings, x, alpha, sparse = sparse)
    penalty_loadings <- checked_pls$loadings
    restore_coef_length <- checked_pls$restore_fun
    x <- checked_pls$trimmed_x
  } else {
    restore_coef_length <- function (coef) coef
    penalty_loadings <- NULL
  }

  if (ncol(x) == 0L) {
    warn("All values in `penalty_loadings` are infinite. Only computing the intercept.")
    std_data <- .standardize_data(matrix(runif(x_dim[[1L]]), ncol = 1L), y, intercept = intercept, sparse = sparse,
                                  standardize = standardize, robust = FALSE)
    lambda <- .elnet_max_lambda(std_data$x, std_data$y, alpha, weights, NULL)

    return(list(std_data = std_data, alpha = alpha, lambda = lambda, weights = weights,
                penalty_loadings = NULL, intercept = intercept, optional_args = optional_args,
                restore_coef_length = restore_coef_length))
  }

  std_data <- .standardize_data(x, y, intercept = intercept, standardize = isTRUE(standardize), robust = FALSE,
                                sparse = sparse)

  # Scale penalty loadings appropriately
  penalty_loadings <- penalty_loadings / std_data$scale_x
  if (length(penalty_loadings) == 0L) {
    penalty_loadings <- NULL
  }

  lambda <- if (missing(lambda) || is.null(lambda)) {
    # Generate lambda sequence automatically.
    if (missing(lambda_min_ratio) || is.null(lambda_min_ratio)) {
      lambda_min_ratio <- max(0.01, alpha) * if (x_dim[[1L]] > x_dim[[2L]]) { 1e-3 } else { 1e-2 }
    }
    max_lambda <- .elnet_max_lambda(std_data$x, std_data$y, alpha, weights, penalty_loadings)
    rev(exp(seq(log(lambda_min_ratio * max_lambda), log(max_lambda), length.out = nlambda)))
  } else {
    sort(.as(lambda, 'numeric'), decreasing = TRUE)
  }

  return(list(std_data = std_data, binary_response = response$binary, alpha = alpha, lambda = lambda,
              weights = weights, penalty_loadings = penalty_loadings, intercept = intercept,
              optional_args = optional_args, restore_coef_length = restore_coef_length))
}

#' @importFrom stats cov
.elnet_max_lambda <- function (x, y, alpha, weights, penalty_loadings) {
  max_cov <- if (is.null(weights)) {
    if (is.null(penalty_loadings)) {
      max(abs(cov(x, y)))
    } else {
      max(abs(cov(x, y)) / penalty_loadings)
    }
  } else {
    if (is.null(penalty_loadings)) {
      max(abs(cov(x, y * weights)))
    } else {
      max(abs(cov(x, y * weights)) / penalty_loadings)
    }
  }

  ((nrow(x) - 1) / nrow(x)) * max_cov / max(0.01, alpha)
}

