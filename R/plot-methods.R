#' Plot Method for Penalized Estimates
#'
#' Plot the coefficient path for fitted penalized elastic net S- or LS-estimates of regression.
#'
#' @param x fitted estimates.
#' @param alpha Plot the coefficient path for the fit with the given hyper-parameter value.
#'   If missing of `NULL`, the first value in `x$alpha` is used.
#' @param ... currently ignored.
#'
#' @family functions for plotting and printing
#'
#' @example examples/pense_fit.R
#' @export
#' @importFrom rlang abort
plot.pense_fit <- function (x, alpha, ...) {
  if (missing(alpha) || is.null(alpha)) {
    alpha <- x$alpha[[1L]]
  }
  alpha <- .as(alpha[[1L]], 'numeric')
  ai <- which((x$alpha - alpha)^2 < .Machine$double.eps)
  if (length(ai) != 1L) {
    abort("Requested `alpha` is not available in the fit.")
  }
  .plot_coef_path(x, lambda_seq = x$lambda[[ai]], alpha = alpha, envir = parent.frame())
  invisible(x)
}

#' Plot Method for Penalized Estimates With Cross-Validation
#'
#' Plot the cross-validation performance or the coefficient path for fitted penalized
#' elastic net S- or LS-estimates of regression.
#'
#' @param x fitted estimates with cross-validation information.
#' @param what plot either the CV performance or the coefficient path.
#' @param alpha If `what = "cv"`, only CV performance for fits with matching `alpha` are plotted.
#'   In case `alpha` is missing or `NULL`, all fits in `x` are plotted.
#'   If `what = "coef.path"`, plot the coefficient path for the fit with the given
#'   hyper-parameter value or, in case `alpha` is missing, for the first value in `x$alpha`.
#' @param se_mult if plotting CV performance, multiplier of the estimated SE.
#' @param ... currently ignored.
#'
#' @family functions for plotting and printing
#' @example examples/pense_fit.R
#' @export
#' @importFrom rlang abort warn
plot.pense_cvfit <- function(x, what = c('cv', 'coef.path'), alpha = NULL, se_mult = 1, ...) {
  what <- match.arg(what)

  if (!any(x$cvres$cvse > 0)) {
    if (identical(what, 'cv') && !missing(se_mult) && isTRUE(se_mult > 0)) {
      warn(paste("Only a single cross-validation replication was performed.",
                 "Standard errors not available."))
    }
    se_mult <- 0
  }

  if (identical(what, 'coef.path')) {
    if (isFALSE(x$call$fit_all)) {
      warn("`x` was created with `fit_all = FALSE`. Coefficient path not available.")
    }
    alpha <- if (is.null(alpha)) {
      x$alpha[[1L]]
    } else {
      .as(alpha[[1L]], 'numeric')
    }
    cvres_rows <- which((x$cvres$alpha - alpha)^2 < .Machine$double.eps)
    if (length(cvres_rows) == 0L) {
      abort("Requested `alpha` not available in the fit `object`.")
    }
    alpha_index <- which.min((x$alpha - alpha)^2)
    se_sel <- .cv_se_selection(x$cvres$cvavg[cvres_rows], x$cvres$cvse[cvres_rows], se_mult)
    lambda_opt <- x$cvres$lambda[cvres_rows][[which(se_sel == 'se_fact')]]

    .plot_coef_path(x, lambda_seq = x$lambda[[alpha_index]],
                    alpha = alpha, lambda_opt = lambda_opt,
                    envir = parent.frame())
  } else {
    .plot_cv_res(x, se_mult = se_mult, alpha = alpha)
  }

  invisible(x)
}

#' @importFrom graphics plot segments abline
#' @importFrom rlang abort warn
.plot_cv_res <- function (object, alpha = NULL, se_mult) {
  measure_label <- switch(object$cv_measure, mape = "Median absolute prediction error",
                          rmspe = "Root mean square prediction error",
                          auroc = "1 - AUROC",
                          tau_size = expression(paste(tau, "-scale of the prediction error")),
                          "Prediction error")
  if (is.null(alpha)) {
    alpha_seq <- object$alpha
  } else {
    alpha_seq <- object$alpha[na.omit(.approx_match(.as(alpha, 'numeric'), object$alpha))]
    if (length(alpha_seq) == 0L) {
      abort("None of the requested `alpha` values is available in the fit")
    }
  }

  colors <- c(none = 'gray20', min = '#0072B2', se_fact = '#56B4E9')
  sizes <- c(none = 1, min = 1.5, se_fact = 1.5)
  shapes <- rep_len(c(16, 15, 17, 18, 0:14), length.out = length(alpha_seq))

  with(object$cvres, {
    ymin <- cvavg - se_mult * cvse
    y_limit_min <- if (any(ymin <= 0)) {
      warn("Error bars extending into the negative range are truncated.")
      if (any(ymin > 0)) {
        min(ymin[ymin > 0])
      } else {
        # All errorbars extend to the negative. Use arbitrary lower limit.
        0.95 * min(cvavg)
      }
    } else {
      min(ymin)
    }
    plot(1, 1, type = 'n', log = 'xy',
         ylim = c(y_limit_min, max(cvavg + se_mult * cvse)),
         xlim = range(lambda),
         xlab = expression(paste('Penalization level (', lambda, ')')),
         ylab = measure_label)

    if (length(alpha_seq) > 1L) {
      title(main = "CV prediction performance")
      legend('topleft', title = expression(alpha), legend = sprintf('%g', alpha_seq),
             pch = shapes, cex = 1)
    } else {
      title(main = sprintf("CV prediction performance for alpha = %g", alpha_seq))
    }

    for (ai in seq_along(alpha_seq)) {
      rows <- which((alpha - alpha_seq[[ai]])^2 < .Machine$double.eps)
      se_sel <- .cv_se_selection(cvm = cvavg[rows], cvsd = cvse[rows], se_fact = se_mult)
      min_ind <- which(se_sel == 'min')
      if (length(min_ind) == 0L) {
        # minimum and SE rule coincide
        min_ind <- which(se_sel == 'se_fact')
      }
      cols <- colors[as.character(se_sel)]
      cex <- sizes[as.character(se_sel)]

      points(x = lambda[rows], y = cvavg[rows], col = cols, cex = cex, pch = shapes[[ai]])

      if (isTRUE(se_mult > 0)) {
        # Ensure errorbars don't extend into the negative!
        errorbar_ymin <- cvavg[rows] - se_mult * cvse[rows]
        neg_errorbars <- which(errorbar_ymin <= 0)
        if (length(neg_errorbars) > 0L) {
          if (length(neg_errorbars) < length(errorbar_ymin)) {
            # There are errorbars that do don't extend to the negative.
            # Use the smallest positive value.
            errorbar_ymin[neg_errorbars] <- min(errorbar_ymin[-neg_errorbars])
          } else {
            # All errorbars extend to the negative. Use arbitrary value.
            errorbar_ymin <- rep.int(cvavg[rows] * 0.95, length(errorbar_ymin))
          }
        }
        segments(x0 = lambda[rows], y0 = errorbar_ymin,
                 x1 = lambda[rows], y1 = cvavg[rows] + se_mult * cvse[rows],
                 col = cols)
      }
    }
  })
}

#' @importFrom graphics plot lines text abline
#' @importFrom rlang abort
.plot_coef_path <- function(object, lambda_seq, alpha, lambda_opt, envir) {
  if (missing(alpha) || is.null(alpha)) {
    alpha <- object$alpha[[1L]]
  }
  alpha <- .as(alpha[[1L]], 'numeric')

  var_names <- tryCatch(colnames(eval(object$call$x, envir = envir)), error = function (e) NULL)
  if (is.null(var_names)) {
    var_names <- paste('X', seq_len(length(object$estimates[[1]]$beta)), sep = '')
  }

  match_ests <- which(vapply(object$estimates, FUN.VALUE = logical(1L), FUN = function (est) {
    (est$alpha - alpha)^2 < .Machine$double.eps
  }))
  if (length(match_ests) == 0L) {
    abort("Requested `alpha` not available in the fit.")
  }

  active_vars <- do.call(rbind, lapply(object$estimates[match_ests], function (est) {
    actives <- .active_indices_and_values(est)
    if (length(actives$var_index) > 0L) {
      return(data.frame(actives, lambda = est$lambda, alpha = est$alpha))
    } else {
      return(data.frame(var_index = integer(), value = numeric(),
                        lambda = numeric(), alpha = numeric()))
    }
  }))

  active_vars$var <- factor(var_names[active_vars$var_index],
                            levels = var_names[sort(unique(active_vars$var_index))])

  var_colors <- (seq_len(nlevels(active_vars$var)) - 1L) %% 10L + 1L
  names(var_colors) <- levels(active_vars$var)

  lambda_seq <- rev(lambda_seq)
  active_vars <- active_vars[order(active_vars$lambda), ]

  plot(1, 0, xlim = range(lambda_seq), ylim = range(0, active_vars$value), type = "n", log = "x",
       xlab = expression(lambda), ylab = "Coefficient estimate",
       main = sprintf("Coefficient path for alpha = %g", alpha))

  for (i in seq_len(nlevels(active_vars$var))) {
    sel_var <- levels(active_vars$var)[[i]]
    with(active_vars, {
      act_values <- value[var == sel_var]
      act_lambda <- lambda[var == sel_var]
      matching_lambdas <- match(act_lambda, lambda_seq)
      lambda_bound <- range(matching_lambdas, na.rm = TRUE)
      pad_0s <- 0L
      if (lambda_bound[[1L]] > 1L) {
        act_lambda <- c(act_lambda, lambda_seq[[lambda_bound[[1L]] - 1L]], min(lambda))
        pad_0s <- 2L
      }
      if (lambda_bound[[2L]] < length(lambda_seq)) {
        act_lambda <- c(act_lambda, lambda_seq[[lambda_bound[[2L]] + 1L]], max(lambda_seq))
        pad_0s <- pad_0s + 2L
      }
      lines(x = act_lambda, y = c(act_values, numeric(pad_0s)), col = var_colors[[sel_var]])
      text(x = min(act_lambda), y = act_values[which.min(act_lambda)], labels = sel_var,
           adj = 1, pos = 2, cex = 0.9, col = "#4f4f4f")
    })
  }

  if (!missing(lambda_opt)) {
    abline(v = lambda_opt, lty = 2, col = "#4f4f4f")
  }
}

#' @importFrom methods is
.active_indices_and_values <- function (est) {
  if (is(est$beta, 'dsparseVector')) {
    return(list(var_index = est$beta@i, value = est$beta@x))
  }
  indices <- which(abs(est$beta) > .Machine$double.eps)
  return(list(var_index = indices, value = est$beta[indices]))
}
