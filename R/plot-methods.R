#' Plot Method for Penalized Estimates
#'
#' Plot the coefficient path for fitted penalized elastic net S- or LS-estimates of regression.
#'
#' @param x fitted estimates.
#' @param ... currently ignored.
#'
#' @family functions for plotting and printing
#'
#' @example examples/pense_fit.R
#' @export
plot.pense_fit <- function (x, ...) {
  .plot_coef_path(x, x$lambda, envir = parent.frame())
}

#' Plot Method for Penalized Estimates With Cross-Validation
#'
#' Plot the cross-validation performance or the coefficient path for fitted penalized elastic net S- or LS-estimates
#' of regression.
#'
#' @param x fitted estimates with cross-validation information.
#' @param what plot either the CV performance or the coefficient path.
#' @param se_mult if plotting CV performance, multiplier of the estimated SE.
#' @param ... currently ignored.
#'
#' @family functions for plotting and printing
#' @example examples/pense_fit.R
#' @export
plot.pense_cvfit <- function(x, what = c('cv', 'coef.path'), se_mult = 1, ...) {
  what <- match.arg(what)
  if (what == 'coef.path' && isFALSE(x$call$fit_all)) {
    stop("`x` was created with `fit_all = FALSE`. Coefficient path not available.")
  } else if (!isTRUE(ncol(x$cv_replications) > 1L)) {
    if (what == 'cv' && !missing(se_mult) && isTRUE(se_mult > 0)) {
      warning("Only a single cross-validation replication was performed. Standard error not available.")
    }
    se_mult <- 0
  }
  se_sel <- .cv_se_selection(x$cvres$cvavg, x$cvres$cvse, se_mult)
  lambda_opt <- x$cvres$lambda[[which(se_sel == 'se_fact')]]
  switch(what, coef.path = .plot_coef_path(x, x$cvres$lambda, lambda_opt, envir = parent.frame()),
         cv = .plot_cv_res(x, se_mult, se_sel))
}

#' @importFrom graphics plot segments abline
#' @importFrom rlang warn
.plot_cv_res <- function (object, se_mult, se_sel) {
  measure_label <- switch(object$cv_measure, mape = "Median absolute prediction error",
                          rmspe = "Root mean square prediction error",
                          auroc = "1 - AUROC",
                          tau_size = expression(paste(tau, "-scale of the prediction error")),
                          "Prediction error")
  colors <- c('gray20', '#0072B2', '#56B4E9')

  with(object$cvres, {
    min_ind <- which(se_sel == 'min')
    if (length(min_ind) == 0L) {
      # minimum and SE rule coincide
      min_ind <- which(se_sel == 'se_fact')
    }
    cols <- colors[as.integer(se_sel)]
    xrange <- range(lambda)

    # Ensure errorbars don't extend into the negative!
    errorbar_ymin <- cvavg - se_mult * cvse
    neg_errorbars <- which(errorbar_ymin <= 0)
    if (length(neg_errorbars) > 0L) {
      warn("Error bars extending into the negative range are truncated.")
      if (length(neg_errorbars) < length(cvavg)) {
        # There are errorbars that do don't extend to the negative. Use the smallest positive value.
        errorbar_ymin[neg_errorbars] <- min(errorbar_ymin[-neg_errorbars])
      } else {
        # All errorbars extend to the negative. Use arbitrary value.
        errorbar_ymin <- rep.int(cvavg * 0.95, length(errorbar_ymin))
      }
    }

    plot(lambda, cvavg, log = 'xy', col = cols,
         ylim = range(cvavg + se_mult * cvse, errorbar_ymin),
         pch = 20L, xlab = expression(paste('Penalization level (', lambda, ')')),
                                      ylab = measure_label, main = "CV prediction performance")
    if (isTRUE(se_mult > 0)) {
      segments(lambda, errorbar_ymin, lambda, cvavg + se_mult * cvse, col = cols)
      abline(h = cvavg[[min_ind]] + se_mult * cvse[[min_ind]], col = colors[[3L]], lty = '22')
    }
  })
}

#' @importFrom graphics plot lines text abline
.plot_coef_path <- function(object, lambda_seq, lambda_opt, envir) {
  var_names <- tryCatch(colnames(eval(object$call$x, envir = envir)), error = function (e) NULL)
  if (is.null(var_names)) {
    var_names <- paste('X', seq_len(length(object$estimates[[1]]$beta)), sep = '')
  }

  active_vars <- do.call(rbind, lapply(object$estimates, function (est) {
    actives <- .active_indices_and_values(est)
    if (length(actives$var_index) > 0L) {
      return(data.frame(.active_indices_and_values(est), lambda = est$lambda))
    } else {
      return(NULL)
    }
  }))
  active_vars$var <- factor(var_names[active_vars$var_index],
                            levels = var_names[sort(unique(active_vars$var_index))])

  var_colors <- (seq_len(nlevels(active_vars$var)) - 1L) %% 10L + 1L
  names(var_colors) <- levels(active_vars$var)

  lambda_seq <- rev(lambda_seq)
  active_vars <- active_vars[order(active_vars$lambda), ]

  plot(1, 0, xlim = range(lambda_seq), ylim = range(0, active_vars$value), type = "n", log = "x",
       xlab = expression(lambda), ylab = "Coefficient estimate")

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
