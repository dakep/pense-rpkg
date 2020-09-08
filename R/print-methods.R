#' Summarize Cross-Validated PENSE Fit
#'
#' If `lambda = "se"` and `object` contains fitted estimates for every penalization level in the sequence, extract the
#' coefficients of the most parsimonious model with prediction performance statistically indistinguishable from the best
#' model. This is determined to be the model with prediction performance within `se_mult * cv_se` from the best model.
#'
#' @param x an (adaptive) PENSE fit with cross-validation information.
#' @param lambda either a string specifying which penalty level to use (`"min"` or `"se"`) or a a single numeric
#'    value of the penalty parameter. See details.
#' @param se_mult If `lambda = "se"`, the multiple of standard errors to tolerate.
#'
#' @seealso [prediction_performance()] for information about the estimated prediction performance.
#' @seealso [coef.cv_pensefit()] for extracting only the estimated coefficients.
#'
#' @importFrom methods is
#' @export
#' @method summary pense_cvfit
summary.pense_cvfit <- function (x, lambda = c('min', 'se'), se_mult = 1, ...) {
  coef_est <- eval.parent(coef(x, lambda, se_mult, sparse = FALSE, add_lambda = TRUE))
  method_name <- if (is(x, 'adapense')) {
    "Adaptive PENSE"
  } else if (is(x, 'pense')) {
    "PENSE"
  } else if (is(x, 'pense_en')) {
    "EN"
  } else if (is(x, 'mest')) {
    "Regularized M"
  }
  tryCatch({
    cat(method_name, "fit", "with", "prediction", "performance", "estimated", "by", x$call$cv_repl, "replications",
        "of", sprintf("%d-fold", x$call$cv_k), "cross-validation.\n", sep = " ", fill = TRUE)
  }, error = function(...) {
    cat(method_name, "fit", "with", "prediction", "performance", "estimated", "via", "cross-validation.\n",
        sep = " ", fill = TRUE)
  })

  nnz_ind <- 1L + which(abs(coef_est[-1L]) > .Machine$double.eps)
  cat(length(nnz_ind), "out", "of", length(coef_est) - 1L, "predictors", "have", "non-zero", "coefficients:\n",
      sep = " ", fill = TRUE)
  print(cbind(Estimate = coef_est[c(1L, nnz_ind)]))
  cat("---\n\n")

  if (is.character(lambda)) {
    lambda_ind <- .lambda_index_cvfit(x, lambda, se_mult)
    lambda <- x$lambda[[lambda_ind]]
    pred_perf <- x$cvres$cvavg[[lambda_ind]]
  } else {
    pred_perf <- NULL
  }

  cat(sprintf("Hyper-parameters: lambda=%s", format(lambda)))
  if (!is.null(x$call$alpha)) {
    cat(sprintf(", alpha=%s", format(x$call$alpha)))
  }
  if (!is.null(x$exponent)) {
    cat(sprintf(", exponent=%s", format(x$exponent)))
  }
  cat("\n")
  if (!is.null(pred_perf)) {
    cat("Estimated", "scale", "of", "the", "prediction", "error:", sprintf("%s\n\n", format(pred_perf)), sep = " ",
        fill = TRUE)
  }
}

#' @rdname summary.pense_cvfit
#' @keywords internal
#' @export
#' @method print pense_cvfit
print.pense_cvfit <- function(x, lambda = c('min', 'se'), se_mult = 1, ...) {
  cl <- match.call(expand.dots = FALSE)
  cl[[1L]] <- quote(summary)
  eval.parent(cl)
}


#' Prediction Performance of Adaptive PENSE Fits
#'
#' Extract the prediction performance of one or more (adaptive) PENSE fits.
#'
#' If `lambda = "se"` and the cross-validation was performed with multiple replications, use the penalty level whit
#' prediction performance within `se_mult` of the best prediction performance.
#'
#' @param ... one or more (adaptive) PENSE fits with cross-validation information.
#' @param lambda a string specifying which penalty level to use (`"min"` or `"se"`) . See details.
#' @param se_mult If `lambda = "se"`, the multiple of standard errors to tolerate.
#'
#' @return a data frame with details about the prediction performance of the given PENSE fits. The data frame
#'    has a custom print method summarizing the prediction performances.
#'
#' @seealso [summarize.pense_cvfit()] for a summary of the fitted model.
#'
#' @importFrom methods is
#' @importFrom rlang abort
#' @export
prediction_performance <- function (..., lambda = c('min', 'se'), se_mult = 1) {
  se_mult <- if (match.arg(lambda) == 'se') {
    .as(se_mult[[1L]], 'numeric')
  } else {
    0
  }

  if (se_mult < 0) {
    abort("`se_mult` must be positive.")
  }

  objects <- as.list(match.call(expand.dots = FALSE)$...)
  object_fallback_names <- sapply(objects, function (obj) paste(deparse(obj), collapse = ' '))
  object_names <- names(objects)
  if (is.null(object_names)) {
    object_names <- object_fallback_names
  } else {
    object_names[object_names == ''] <- object_fallback_names[object_names == '']
  }
  names(objects) <- object_names

  pred_perf <- do.call(rbind, lapply(object_names, function (on) {
    object <- eval.parent(objects[[on]], n = 2)
    if (!is(object, 'pense_cvfit')) {
      abort(sprintf("`%s` must be a cross-validated PENSE fit.", on))
    }
    if (!isTRUE(ncol(object$cv_replications) > 1L) && isTRUE(se_mult > 0)) {
      warn(sprintf(
        "Only a single cross-validation replication was performed for object `%s`. Standard error not available."), on)
      se_mult <- 0
    }
    se_selection <- which(.cv_se_selection(object$cvres$cvavg, object$cvres$cvse, se_mult) == 'se_fact')
    cvres <- object$cvres[se_selection, ]

    cvres$model_size <- if (isFALSE(object$call$fit_all)) {
      if (se_mult > 0) {
        NA_integer_
      } else {
        sum(abs(object$estimates[[1L]]$beta) > .Machine$double.eps)
      }
    } else {
      sum(abs(object$estimates[[se_selection]]$beta) > .Machine$double.eps)
    }
    cvres$name <- on
    cvres$alpha <- if (!is.null(object$call$alpha)) { object$call$alpha } else { NA_real_ }
    cvres$exponent <- if (!is.null(object$exponent)) { object$exponent } else { NA_real_ }

    return(cvres)
  }))
  pred_perf <- pred_perf[order(pred_perf$cvavg), c('name', 'cvavg', 'cvse', 'model_size', 'lambda', 'alpha',
                                                   'exponent')]
  rownames(pred_perf) <- NULL
  class(pred_perf) <- c('pense_pred_perf', class(pred_perf))
  return(pred_perf)
}

#' @rdname prediction_performance
#' @export
#' @method print pense_pred_perf
print.pense_pred_perf <- function (x, ...) {
  xprint <- x[ , c('name', 'cvavg', 'cvse', 'model_size', 'alpha', 'exponent')]

  class(xprint) <- 'data.frame'
  colnames(xprint) <- c('Model', 'Estimate', 'Std. Error', 'Predictors', 'alpha', 'exp.')

  if (all(is.na(x$cvse))) {
    xprint$S.E. <- NULL
  }
  if (all(is.na(x$alpha))) {
    xprint$alpha <- NULL
  }
  if (all(is.na(x$exponent))) {
    xprint$exponent <- NULL
  }
  cat("Prediction performance estimated by cross-validation:\n\n")
  print(xprint)
  return(invisible(x))
}

#' Print Metrics
#'
#' Pretty-print a list of metrics from optimization algorithm (if `pense` was built with metrics enabled).
#'
#' @param x metrics object for printing.
#' @param max_level maximum level of printing which is applied for printing nested metrics.
#' @keywords internal
#' @export
#' @method print nsoptim_metrics
print.nsoptim_metrics <- function (x, max_level = NA, ...) {
  .print_metrics(x, max_level, '')
  invisible(NULL)
}

.print_metrics <- function (metrics, max_level, prefix) {
  cat(prefix, '* ', metrics$name, sep = '')

  other_metrics <- setdiff(names(metrics), c('name', 'sub_metrics'))
  if (length(other_metrics) > 0L) {
    cat(':', sep = '')
  }

  for (metric_name in other_metrics) {
    if (is.numeric(metrics[[metric_name]])) {
      cat(sprintf(' %s=%g;', metric_name, metrics[[metric_name]]))
    } else if (is.character(metrics[[metric_name]])) {
      cat(sprintf(' %s="%s";', metric_name,
                  sub('(\\s|;)+$', '', metrics[[metric_name]])))
    } else {
      cat(sprintf(' %s=%s;', metric_name, metrics[[metric_name]]))
    }
  }
  cat('\n', sep = '')

  if (!isTRUE(max_level <= 0L) && !is.null(metrics$sub_metrics)) {
    lapply(rev(metrics$sub_metrics), .print_metrics, max_level = max_level - 1L,
           prefix = paste0(prefix, '  '))
  }
  invisible(NULL)
}
