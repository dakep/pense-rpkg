#' Extract Coefficient Estimates
#'
#' Extract coefficients from a PENSE (or LS-EN) regularization path fitted by [pense()] or [elnet()].
#'
#' @param object PENSE regularization path to extract coefficients from.
#' @param lambda a single value of the penalty parameter.
#' @param sparse should coefficients be returned as sparse or dense vectors? Defaults to the `sparse` argument
#'    supplied to [pense()]. Can also be set to `sparse = 'matrix'`, in which case a sparse matrix
#'    is returned instead of a sparse vector.
#' @param exact defunct Always gives a warning if `lambda` is not part of the fitted sequence and hence coefficients
#'    are interpolated.
#' @param correction defunct.
#' @param ... currently not used.
#' @return either a numeric vector or a sparse vector of type [dsparseVector][Matrix::sparseVector-class]
#'    of size \eqn{p + 1}, depending on the `sparse` argument.
#'    Note: prior to version 2.0.0 sparse coefficients were returned as sparse matrix of type *dgCMatrix*.
#'    To get a sparse matrix, use `sparse = 'matrix'`.
#'
#' @family functions for extracting components
#' @seealso [coef.pense_cvfit()] for extracting coefficients from a PENSE fit with hyper-parameters chosen by
#'    cross-validation
#' @example examples/pense_fit.R
#' @export
#' @importFrom lifecycle deprecate_warn deprecated is_present
#' @importFrom rlang warn
coef.pense_fit <- function (object, lambda, sparse = NULL, exact = deprecated(), correction = deprecated(), ...) {
  if (is_present(exact)) {
    deprecate_warn('2.0.0', 'coef(exact=)')
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
    return(.concat_coefs(object$estimates[[1L]], object$call, sparse, parent.frame(), concat))
  }

  lambda_match <- .approx_match(lambda, object$lambda)
  if (is.na(lambda_match)) {
    warn("Requested penalization level not part of the sequence. Returning interpolated coefficients.")
    return(.interpolate_coefs(object, lambda, object$lambda, sparse = sparse, envir = parent.frame(), concat = concat))
  } else {
    return(.concat_coefs(object$estimates[[lambda_match]], object$call, sparse, parent.frame(), concat))
  }
}


#' Extract Coefficient Estimates
#'
#' Extract coefficients from a PENSE (or LS-EN) regularization path with hyper-parameters chosen by cross-validation.
#'
#' If `lambda = "se"` and `object` contains fitted estimates for every penalization level in the sequence, extract the
#' coefficients of the most parsimonious model with prediction performance statistically indistinguishable from the best
#' model. This is determined to be the model with prediction performance within `se_mult * cv_se` from the best model.
#'
#' @param object PENSE with cross-validated hyper-parameters to extract coefficients from.
#' @param lambda either a string specifying which penalty level to use (`"min"` or `"se"`) or a a single numeric
#'    value of the penalty parameter. See details.
#' @param se_mult If `lambda = "se"`, the multiple of standard errors to tolerate.
#' @param sparse should coefficients be returned as sparse or dense vectors? Defaults to the `sparse` argument
#'    supplied to [pense_cv()]. Can also be set to `sparse = 'matrix'`, in which case a sparse matrix
#'    is returned instead of a sparse vector.
#' @param exact deprecated. Always gives a warning if `lambda` is not part of the fitted sequence and coefficients
#'    are interpolated.
#' @param correction defunct.
#' @param ... currently not used.
#' @return either a numeric vector or a sparse vector of type [dsparseVector][Matrix::sparseVector-class]
#'    of size \eqn{p + 1}, depending on the `sparse` argument.
#'    Note: prior to version 2.0.0 sparse coefficients were returned as sparse matrix of type *dgCMatrix*.
#'    To get a sparse matrix, use `sparse = 'matrix'`.
#'
#' @family functions for extracting components
#' @example examples/pense_fit.R
#' @export
coef.pense_cvfit <- function (object, lambda = c('min', 'se'), se_mult = 1, sparse = NULL, exact = deprecated(),
                              correction = deprecated(), ...) {
  if (is_present(exact)) {
    deprecate_warn('2.0.0', 'coef(exact=)')
  }
  if (is_present(correction)) {
    deprecate_stop('2.0.0', 'coef(correction=)')
  }

  cl <- match.call(expand.dots = TRUE)
  concat <- !isFALSE(cl$concat)
  lambda_index <- .lambda_index_cvfit(object, lambda, se_mult)

  if (is.na(lambda_index)) {
    warn("Requested penalization level not part of the sequence. Returning interpolated coefficients.")
    .interpolate_coefs(object, lambda, object$lambda, sparse = sparse, envir = parent.frame(), concat = concat)
  } else {
    .concat_coefs(object$estimates[[lambda_index]], object$call, sparse, parent.frame(), concat)
  }
}

#' @importFrom rlang warn
.lambda_index_cvfit <- function (object, lambda = c('min', 'se'), se_mult) {
  if (is.character(lambda)) {
    lambda <- match.arg(lambda)
  }

  if (isFALSE(object$call$fit_all)) {
    if (!isTRUE(lambda == 'min')) {
      warn(paste("`object` was created with `fit_all = FALSE`. Only the estimate at the minimum is available and will",
                 "be returned."))
    }
    return(1L)
  }

  lambda <- if (is.character(lambda)) {
    if (!isTRUE(ncol(object$cv_replications) > 1L) && isTRUE(lambda == 'se') && isTRUE(se_mult > 0)) {
      warn(paste("Only a single cross-validation replication was performed. Standard error not available.",
                 "Using minimum lambda."))
    }

    if (isTRUE(lambda == 'min')) {
      lambda <- 'se'
      se_mult <- 0
    }

    se_selection <- .cv_se_selection(object$cvres$cvavg, object$cvres$cvse, se_mult)
    return(which(se_selection == 'se_fact'))
  } else {
    if (length(lambda) > 1L) {
      warn("Only first element in `lambda` is used.")
    }
    .as(lambda[[1L]], 'numeric')
  }

  if (isTRUE(lambda > object$lambda[[1L]]) && isTRUE(sum(abs(object$estimates[[1L]]$beta)) < .Machine$double.eps)) {
    return(1L)
  }

  return(.approx_match(lambda, object$lambda))
}

## Interpolate the coefficients at `lambda` using estimates from `object` at `lambda_seq`.
## @param ... passed on to `.concat_coefs()`
.interpolate_coefs <- function (object, lambda, lambda_seq, ...) {
  if (lambda < lambda_seq[[length(lambda_seq)]]) {
    return(.concat_coefs(object$estimates[[length(lambda_seq)]], object$call))
  } else if (lambda > lambda_seq[[1L]]) {
    return(.concat_coefs(object$estimates[[1L]], object$call))
  }

  lambda_diff <- lambda_seq - lambda
  lambda_diff_abs <- abs(lambda_diff)
  # take the average of the closest two solutions
  lambda_ind_left <- min(which(lambda_diff < 0))
  lambda_ind_right <- max(which(lambda_diff > 0))
  interp_w <- 1 / lambda_diff_abs[c(lambda_ind_left, lambda_ind_right)]
  interp_beta <- (interp_w[[1L]] * object$estimates[[lambda_ind_left]]$beta +
                    interp_w[[2L]] * object$estimates[[lambda_ind_right]]$beta) / sum(interp_w)
  interp_int <- (interp_w[[1L]] * object$estimates[[lambda_ind_left]]$intercept +
                   interp_w[[2L]] * object$estimates[[lambda_ind_right]]$intercept) / sum(interp_w)
  return(.concat_coefs(list(intercept = interp_int, beta = interp_beta), object$call, ...))
}

#' @importFrom methods is
#' @importFrom Matrix sparseVector sparseMatrix
.concat_coefs <- function (coef, call, sparse, envir, concat = TRUE) {
  if (!concat) {
    return(coef)
  }

  # Determine names
  var_names <- tryCatch(colnames(eval(call$x, envir = envir)), error = function (e) NULL)
  if (is.null(var_names)) {
    var_names <- paste('X', seq_len(length(coef$beta)), sep = '')
  }
  var_names <- c('(Intercept)', var_names)

  if (!is.null(sparse) && pmatch(sparse[[1L]], 'matrix', nomatch = 0L) == 1L) {
    if (is(coef$beta, 'dsparseVector')) {
      return(sparseMatrix(x = c(coef$intercept, coef$beta@x), i = c(1L, 1L + coef$beta@i),
                          j = rep.int(1L, length(coef$beta@i) + 1L),
                          dims = c(coef$beta@length + 1L, 1L), dimnames = list(var_names, NULL)))
    } else {
      nz_ind <- which(abs(coef$beta) > .Machine$double.eps)
      return(sparseMatrix(x = c(coef$intercept, coef$beta[nz_ind]), i = c(1L, 1L + nz_ind),
                          j = rep.int(1L, length(nz_ind) + 1L),
                          dims = c(length(coef$beta) + 1L, 1L), dimnames = list(var_names, NULL)))
    }
  } else if (isTRUE(sparse) || (is.null(sparse) && is(coef$beta, 'dsparseVector'))) {
    if (is(coef$beta, 'dsparseVector')) {
      return(sparseVector(c(coef$intercept, coef$beta@x), i = c(1L, coef$beta@i + 1L), length = coef$beta@length + 1L))
    } else {
      return(sparseVector(c(coef$intercept, coef$beta), i = seq_len(length(coef$beta) + 1L),
                          length = length(coef$beta) + 1L))
    }
  } else if (isFALSE(sparse) || (is.null(sparse) && is.numeric(coef$beta))) {
    coefvec <- c(coef$intercept, as.numeric(coef$beta))
    names(coefvec) <- var_names
    return(coefvec)
  } else {
    abort("argument `sparse` must be TRUE/FALSE or \"matrix\".")
  }
}
