#' Extract Coefficient Estimates
#'
#' Extract coefficients from an adaptive PENSE (or LS-EN) regularization path fitted by [pense()]
#' or [elnet()].
#'
#' @param object PENSE regularization path to extract coefficients from.
#' @param lambda a single value of the penalty parameter.
#' @param alpha Either a single number or missing.
#'    If given, only fits with the given `alpha` value are considered.
#'    If `object` was fit with multiple `alpha` values, the parameter `alpha` must not be missing.
#' @param sparse should coefficients be returned as sparse or dense vectors? Defaults to the
#'    sparsity setting in `object`.
#'    Can also be set to `sparse = 'matrix'`, in which case a sparse matrix
#'    is returned instead of a sparse vector.
#' @param standardized return the standardized coefficients.
#' @param exact,correction defunct.
#' @param ... currently not used.
#' @return either a numeric vector or a sparse vector of type
#'    [dsparseVector][Matrix::sparseVector-class]
#'    of size \eqn{p + 1}, depending on the `sparse` argument.
#'    Note: prior to version 2.0.0 sparse coefficients were returned as sparse matrix
#'    of type *dgCMatrix*.
#'    To get a sparse matrix as in previous versions, use `sparse = 'matrix'`.
#'
#' @family functions for extracting components
#' @seealso [coef.pense_cvfit()] for extracting coefficients from a PENSE fit with
#'    hyper-parameters chosen by cross-validation
#' @example examples/pense_fit.R
#' @export
#' @importFrom lifecycle deprecate_warn deprecated is_present
#' @importFrom rlang warn
coef.pense_fit <- function (object, lambda, alpha, sparse = NULL, standardized = FALSE,
                            exact = deprecated(), correction = deprecated(), ...) {
  if (is_present(exact)) {
    deprecate_stop('2.0.0', 'coef(exact=)')
  }
  if (is_present(correction)) {
    deprecate_stop('2.0.0', 'coef(correction=)')
  }

  if (length(lambda) > 1L) {
    warn("Only first element in `lambda` is used.")
  }

  cl <- match.call(expand.dots = TRUE)
  concat <- !isFALSE(cl$concat)

  lambda <- .as(lambda[[1L]], 'numeric')

  if (isTRUE(lambda > object$lambda[[1L]]) && isTRUE(sum(abs(object$estimates[[1L]]$beta)) < .Machine$double.eps)) {
    return(.concat_coefs(object$estimates[[1L]], object$call, sparse = sparse, envir = parent.frame(),
                         standardized = standardized, concat = concat))
  }

  lambda_index <- .lambda_index_cvfit(object, lambda = lambda, alpha = alpha, se_mult = se_mult)
  if (length(lambda_index) > 1L) {
    warn(paste("Requested penalization level not part of the sequence.",
               "Returning interpolated coefficients."))
    return(.interpolate_coefs(object, indices = lambda_index, lambda = lambda,
                              sparse = sparse, envir = parent.frame(),
                              standardized = standardized, concat = concat))
  } else {
    return(.concat_coefs(object$estimates[[lambda_index]], object$call,
                         sparse = sparse, envir = parent.frame(),
                         standardized = standardized, concat = concat))
  }
}

#' Extract Coefficient Estimates
#'
#' Extract coefficients from an adaptive PENSE (or LS-EN) regularization path with hyper-parameters
#' chosen by cross-validation.
#'
#' If `lambda = "se"` and `object` contains fitted estimates for every penalization level in
#' the sequence, extract the coefficients of the most parsimonious model with prediction
#' performance statistically indistinguishable from the best model.
#' This is determined to be the model with prediction performance within `se_mult * cv_se`
#' from the best model.
#'
#' @param object PENSE with cross-validated hyper-parameters to extract coefficients from.
#' @param lambda either a string specifying which penalty level to use (`"min"` or `"se"`) or
#'    a single numeric value of the penalty parameter. See details.
#' @param alpha Either a single number or missing.
#'    If given, only fits with the given `alpha` value are considered.
#'    If `lambda` is a numeric value and `object` was fit with multiple `alpha`
#'    values, the parameter `alpha` must not be missing.
#' @param se_mult If `lambda = "se"`, the multiple of standard errors to tolerate.
#' @param sparse should coefficients be returned as sparse or dense vectors?
#'    Defaults to the sparsity setting of the given `object`.
#'    Can also be set to `sparse = 'matrix'`, in which case a sparse matrix
#'    is returned instead of a sparse vector.
#' @param standardized return the standardized coefficients.
#' @param exact,correction defunct.
#' @param ... currently not used.
#' @return either a numeric vector or a sparse vector of type
#'    [dsparseVector][Matrix::sparseVector-class]
#'    of size \eqn{p + 1}, depending on the `sparse` argument.
#'    Note: prior to version 2.0.0 sparse coefficients were returned as sparse matrix of
#'    type *dgCMatrix*.
#'    To get a sparse matrix as in previous versions, use `sparse = 'matrix'`.
#'
#' @family functions for extracting components
#' @example examples/pense_fit.R
#' @export
coef.pense_cvfit <- function (object, lambda = c('min', 'se'), alpha, se_mult = 1, sparse = NULL,
                              standardized = FALSE,
                              exact = deprecated(), correction = deprecated(), ...) {
  if (is_present(exact)) {
    deprecate_stop('2.0.0', 'coef(exact=)')
  }
  if (is_present(correction)) {
    deprecate_stop('2.0.0', 'coef(correction=)')
  }

  cl <- match.call(expand.dots = TRUE)
  concat <- !isFALSE(cl$concat)
  lambda_index <- .lambda_index_cvfit(object, lambda = lambda, alpha = alpha, se_mult = se_mult)

  if (length(lambda_index) > 1L) {
    warn(paste("Requested penalization level not part of the sequence.",
               "Returning interpolated coefficients."))
    .interpolate_coefs(object, indices = lambda_index,
                       lambda = .as(lambda[[1L]], 'numeric'),
                       sparse = sparse, envir = parent.frame(),
                       standardized = standardized, concat = concat)
  } else {
    .concat_coefs(object$estimates[[lambda_index]], call = object$call, sparse = sparse,
                  envir = parent.frame(), standardized = standardized, concat = concat)
  }
}

## Determine the appropriate indices in the `object$estimates` list
#' @importFrom rlang warn abort
.lambda_index_cvfit <- function (object, lambda = c('min', 'se'), alpha, se_mult) {
  if (is.character(lambda)) {
    lambda <- match.arg(lambda)
  }

  if (isFALSE(object$call$fit_all)) {
    if (!isTRUE(lambda == 'min')) {
      warn(paste("`object` was created with `fit_all = FALSE`. Only the estimate at the minimum",
                 "is available and will be returned."))
    }
    return(1L)
  }

  if (is.character(lambda)) {
    if (!any(object$cvres$cvse > 0) && isTRUE(lambda == 'se') && isTRUE(se_mult > 0)) {
      warn(paste("Only a single cross-validation replication was performed.",
                 "Standard error not available. Using minimum lambda."))
    }

    if (isTRUE(lambda == 'min')) {
      lambda <- 'se'
      se_mult <- 0
    }

    considered_alpha <- if (!missing(alpha) && !is.null(alpha)) {
      object$alpha[na.omit(.approx_match(.as(alpha, 'numeric'), object$alpha))]
    } else {
      object$alpha
    }

    if (length(considered_alpha) == 0L) {
      abort("`object` was not fit with the requested `alpha` value.")
    }

    best_per_alpha <- vapply(considered_alpha, FUN.VALUE = numeric(2L), FUN = function (alpha) {
      rows <- which((object$cvres$alpha - alpha)^2 < .Machine$double.eps)
      se_selection <- .cv_se_selection(object$cvres$cvavg[rows],
                                       object$cvres$cvse[rows], se_mult)
      sel_index <- rows[[which(se_selection == 'se_fact')]]
      c(lambda = object$cvres$lambda[[sel_index]],
        cvavg = object$cvres$cvavg[[sel_index]])
    })
    best_alpha_index <- which.min(best_per_alpha['cvavg', ])

    alpha <- considered_alpha[[best_alpha_index]]
    lambda <- best_per_alpha['lambda', best_alpha_index]
  }

  if (length(lambda) > 1L) {
    warn("Only first element in `lambda` is used.")
  }

  if (length(alpha) > 1L) {
    warn("Only the first element in `alpha` is used.")
  }

  if (length(object$alpha) > 1L && (missing(alpha) || is.null(alpha))) {
    abort(paste("`object` was fit for multiple `alpha` values.",
                "Use the `alpha` parameter to select one of these values."))
  }

  alpha <- .as(alpha[[1L]], 'numeric')
  est_alpha <- vapply(object$estimates, FUN = `[[`, 'alpha', FUN.VALUE = numeric(1L))
  alpha_ests_indices <- which((alpha - est_alpha)^2 < .Machine$double.eps)
  if (length(alpha_ests_indices) == 0L) {
    abort("`object` was not fit with the requested `alpha` value.")
  }

  lambda <- .as(lambda[[1L]], 'numeric')
  est_lambda <- vapply(object$estimates[alpha_ests_indices],
                       FUN = `[[`, 'lambda', FUN.VALUE = numeric(1L))

  # Determine the closest match in the lambda grid
  selected_ests <- alpha_ests_indices[.approx_match(lambda, est_lambda)]

  if (anyNA(selected_ests)) {
    # No lambda matches exactly. Use the two surrounding values (if available).
    est_lambda_ord <- order(est_lambda, decreasing = TRUE)
    max_lambda_index <- est_lambda_ord[[1L]]
    min_lambda_index <- est_lambda_ord[[length(est_lambda_ord)]]
    if (!isTRUE(lambda < est_lambda[[max_lambda_index]])) {
      # If lambda is greater than the highest level return the corresponding index.
      selected_ests <- alpha_ests_indices[[max_lambda_index]]
      if (!isTRUE(sum(abs(object$estimates[[selected_ests]]$beta)) < .Machine$double.eps)) {
        warn(paste("Selected `lambda` is larger than highest penalization level in `object`.",
                   "Returning estimate for highest penalization level."))
      }
    } else if (!isTRUE(lambda > est_lambda[[min_lambda_index]])) {
      # If lambda is smaller than the smallest level return the corresponding index.
      warn(paste("Selected `lambda` is smaller than smallest penalization level in `object`.",
                 "Returning estimate for smallest penalization level."))
      selected_ests <- alpha_ests_indices[[min_lambda_index]]
    } else {
      lambda_diff <- est_lambda[est_lambda_ord] - lambda
      lambda_ind_left <- min(which(lambda_diff < 0))
      lambda_ind_right <- max(which(lambda_diff >= 0))
      selected_ests <- alpha_ests_indices[est_lambda_ord[c(lambda_ind_left, lambda_ind_right)]]
    }
  }

  selected_ests
}

## Interpolate the coefficients at `lambda` using estimates from `object` at `lambda_seq`.
## @param ... passed on to `.concat_coefs()`
.interpolate_coefs <- function (object, indices, lambda, ...) {
  lambdas <- vapply(object$estimates[indices], FUN = `[[`, 'lambda', FUN.VALUE = numeric(1L))
  interp_w <- 1 / abs(lambdas - lambda)
  interp_w <- interp_w / sum(interp_w)

  interp_beta <- interp_w[[1L]] * object$estimates[[indices[[1L]]]]$beta +
    interp_w[[2L]] * object$estimates[[indices[[2L]]]]$beta
  interp_int <- interp_w[[1L]] * object$estimates[[indices[[1L]]]]$intercept +
    interp_w[[2L]] * object$estimates[[indices[[2L]]]]$intercept

  interp_std_beta <- interp_w[[1L]] * object$estimates[[indices[[1L]]]]$std_beta +
    interp_w[[2L]] * object$estimates[[indices[[2L]]]]$std_beta
  interp_std_int <- interp_w[[1L]] * object$estimates[[indices[[1L]]]]$std_intercept +
    interp_w[[2L]] * object$estimates[[indices[[2L]]]]$std_intercept

  return(.concat_coefs(list(intercept = interp_int, beta = interp_beta,
                            std_intercept = interp_std_int, std_beta = interp_std_beta),
                       object$call, ...))
}

#' @importFrom methods is
#' @importFrom Matrix sparseVector sparseMatrix
.concat_coefs <- function (coef, call, sparse, envir, standardized = FALSE, concat = TRUE) {
  if (!concat) {
    return(coef)
  }

  # Determine names
  var_names <- tryCatch(colnames(eval(call$x, envir = envir)), error = function (e) NULL)

  if (is.null(var_names)) {
    var_names <- paste('X', seq_len(length(coef$beta)), sep = '')
  }

  var_names <- c('(Intercept)', var_names)

  if (isTRUE(standardized)) {
    beta <- coef$std_beta
    intercept <- coef$std_intercept
  } else {
    beta <- coef$beta
    intercept <- coef$intercept
  }

  if (!is.null(sparse) && pmatch(sparse[[1L]], 'matrix', nomatch = 0L) == 1L) {
    if (is(beta, 'dsparseVector')) {
      return(sparseMatrix(x = c(intercept, beta@x), i = c(1L, 1L + beta@i),
                          j = rep.int(1L, length(beta@i) + 1L),
                          dims = c(beta@length + 1L, 1L), dimnames = list(var_names, NULL)))
    } else {
      nz_ind <- which(abs(beta) > .Machine$double.eps)
      return(sparseMatrix(x = c(intercept, beta[nz_ind]), i = c(1L, 1L + nz_ind),
                          j = rep.int(1L, length(nz_ind) + 1L),
                          dims = c(length(beta) + 1L, 1L), dimnames = list(var_names, NULL)))
    }
  } else if (isTRUE(sparse) || (is.null(sparse) && is(beta, 'dsparseVector'))) {
    if (is(beta, 'dsparseVector')) {
      return(sparseVector(c(intercept, beta@x), i = c(1L, beta@i + 1L), length = beta@length + 1L))
    } else {
      return(sparseVector(c(intercept, beta), i = seq_len(length(beta) + 1L),
                          length = length(beta) + 1L))
    }
  } else if (isFALSE(sparse) || (is.null(sparse) && is.numeric(beta))) {
    coefvec <- c(intercept, as.numeric(beta))
    names(coefvec) <- var_names
    return(coefvec)
  } else {
    abort("argument `sparse` must be TRUE/FALSE or \"matrix\".")
  }
}
