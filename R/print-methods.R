#' Summarize Cross-Validated PENSE Fit
#'
#' If `lambda = "se"` and `object` contains fitted estimates for every penalization level in the sequence, extract the
#' coefficients of the most parsimonious model with prediction performance statistically indistinguishable from the best
#' model. This is determined to be the model with prediction performance within `se_mult * cv_se` from the best model.
#'
#' @param object,x an (adaptive) PENSE fit with cross-validation information.
#' @param lambda either a string specifying which penalty level to use
#'    (`"min"`, `"se"`, `"{x}-se`")
#'    or a single numeric value of the penalty parameter. See details.
#' @param alpha Either a single number or missing.
#'    If given, only fits with the given `alpha` value are considered.
#'    If `lambda` is a numeric value and `object` was fit with multiple `alpha`
#'    values, the parameter `alpha` must not be missing.
#' @param se_mult If `lambda = "se"`, the multiple of standard errors to tolerate.
#' @param ... ignored.
#'
#' @family functions for plotting and printing
#' @seealso [prediction_performance()] for information about the estimated prediction performance.
#' @seealso [coef.pense_cvfit()] for extracting only the estimated coefficients.
#'
#' @export
#' @method summary pense_cvfit
#'
#' @importFrom methods is
#' @importFrom stats coef
summary.pense_cvfit <- function (object, alpha, lambda = 'min', se_mult = 1, ...) {
  coef_est <- eval.parent(coef(object, alpha = alpha, lambda = lambda, se_mult = se_mult,
                               concat = FALSE))
  method_name <- if (is(object, 'adapense')) {
    "Adaptive PENSE"
  } else if (is(object, 'pense')) {
    "PENSE"
  } else if (is(object, 'pense_en')) {
    "EN"
  } else if (is(object, 'mest')) {
    "Regularized M"
  }
  tryCatch({
    cat(method_name, "fit", "with", "prediction", "performance", "estimated", "by",
        object$call$cv_repl, "replications", "of", sprintf("%d-fold", object$call$cv_k),
        "cross-validation.\n", sep = " ", fill = TRUE)
  }, error = function(...) {
    cat(method_name, "fit", "with", "prediction", "performance", "estimated", "via",
        "cross-validation.\n", sep = " ", fill = TRUE)
  })

  coef_est_named <- eval.parent(coef(object, alpha = alpha, lambda = lambda, se_mult = se_mult,
                                     sparse = FALSE))
  nnz_ind <- 1L + which(abs(coef_est_named[-1L]) > .Machine$double.eps)
  cat(length(nnz_ind), "out", "of", length(coef_est$beta), "predictors", "have", "non-zero",
      "coefficients:\n", sep = " ", fill = TRUE)
  print(cbind("Estimate" = coef_est_named[c(1L, nnz_ind)]))
  cat("---\n\n")

  # Try to determine the prediction performance
  cv_res_ind <- which((object$cvres$lambda - coef_est$lambda)^2 < .Machine$double.eps &
                        (object$cvres$alpha - coef_est$alpha)^2 < .Machine$double.eps)
  pred_perf <- if (length(cv_res_ind) > 1L) {
    object$cvres$cvavg[[cv_res_ind[[1L]]]]
  } else {
    NULL
  }

  cat("Hyper-parameters: lambda=", format(coef_est$lambda), sep = '')
  if (!is.null(object$alpha)) {
    cat(", alpha=", format(coef_est$alpha), sep = '')
  }
  if (!is.null(object$exponent)) {
    cat(", exponent=", format(object$exponent), sep = '')
  }
  cat("\n")
  if (!is.null(pred_perf)) {
    cat("Estimated", "scale", "of", "the", "prediction", "error:",
        sprintf("%s\n\n", format(pred_perf)), sep = " ", fill = TRUE)
  }
  invisible(object)
}

#' @rdname summary.pense_cvfit
#' @export
#' @method print pense_cvfit
print.pense_cvfit <- function(x, alpha, lambda = 'min', se_mult = 1, ...) {
  cl <- match.call(expand.dots = FALSE)
  cl[[1L]] <- quote(summary)
  names(cl)[[which(names(cl) == 'x')]] <- 'object'
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
#' @param lambda either a string specifying which penalty level to use
#'    (`"min"`, `"se"`, `"{x}-se`")
#'    or a single numeric value of the penalty parameter. See details.
#' @param alpha Either a numeric vector or `NULL` (default).
#'    If given, only fits with the given `alpha` value are considered.
#'    If `lambda` is a numeric value and `object` was fit with multiple `alpha`
#'    values, the parameter `alpha` must not be missing.
#' @param se_mult If `lambda = "se"`, the multiple of standard errors to tolerate.
#'
#' @return a data frame with details about the prediction performance of the given PENSE fits. The data frame
#'    has a custom print method summarizing the prediction performances.
#'
#' @family functions for plotting and printing
#' @seealso [summary.pense_cvfit()] for a summary of the fitted model.
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom rlang abort
prediction_performance <- function (..., alpha = NULL, lambda = 'min', se_mult = 1) {
  se_mult <- if (is.character(lambda)) {
    if (identical(lambda, 'se')) {
      .as(se_mult[[1L]], 'numeric')
    } else {
      .parse_se_string(lambda, only_fact = TRUE)
    }
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
  eval_frame <- parent.frame()

  pred_perf <- do.call(rbind, lapply(object_names, function (on) {
    object <- eval(objects[[on]], eval_frame)
    if (!is(object, 'pense_cvfit')) {
      abort(sprintf("`%s` must be a cross-validated fit.", on))
    }
    if (isTRUE(se_mult > 0) && !any(object$cvres$cvse > 0)) {
      warn(paste("Only a single cross-validation replication was performed for object `", on,
                 "`. Standard error not available.", sep = ""))
      se_mult <- 0
    }

    sel_alpha <- if (is.null(alpha)) {
      object$alpha
    } else {
      object$alpha[na.omit(.approx_match(alpha, object$alpha))]
    }

    sel_indices <- vapply(sel_alpha, FUN.VALUE = integer(1L), FUN = function (alpha) {
      rows <- which((object$cvres$alpha - alpha)^2 < .Machine$double.eps)
      se_sel <- .cv_se_selection(object$cvres$cvavg[rows], object$cvres$cvse[rows], se_mult)
      rows[[which(se_sel == 'se_fact')[[1L]]]]
    })

    cvres <- object$cvres[sel_indices, ]
    cvres$model_size <- vapply(sel_indices, FUN.VALUE = integer(1L), FUN = function (ind) {
      est_ind <- which(vapply(object$estimates, FUN.VALUE = logical(1L), FUN = function (est) {
        (est$alpha - object$cvres$alpha[[ind]])^2 < .Machine$double.eps &
          (est$lambda - object$cvres$lambda[[ind]])^2 < .Machine$double.eps
      }))
      if (length(est_ind) > 0L) {
        sum(abs(object$estimates[[est_ind]]$beta) > .Machine$double.eps)
      } else {
        NA_integer_
      }
    })
    cvres$name <- rep.int(on, length(sel_indices))
    cvres$exponent <- rep.int(if (!is.null(object$exponent)) { object$exponent } else { NA_real_ },
                              length(sel_indices))

    return(cvres)
  }))
  pred_perf <- pred_perf[order(pred_perf$cvavg),
                         c('name', 'cvavg', 'cvse', 'model_size', 'lambda', 'alpha', 'exponent')]
  rownames(pred_perf) <- NULL
  class(pred_perf) <- c('pense_pred_perf', class(pred_perf))
  return(pred_perf)
}

#' @rdname prediction_performance
#' @param x an object with information on prediction performance created with `prediction_performance()`.
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
  invisible(metrics)
}
