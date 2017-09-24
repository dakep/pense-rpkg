#' Plot Method for Fitted Penalized Elastic Net S/MM-Estimates of Regression
#'
#' Plot the cross-validation error or the coefficient path for a fitted
#' PENSE estimate.
#'
#' @param x a PENSE or PENSEM estimate from \code{\link{pense}} or
#'      \code{\link{pensem}}.
#' @param what plot either the cross-validated prediction error
#'      (\code{"cv"}; default) or the coefficient paths.
#' @param ... currently ignored.
#'
#' @example examples/pense-plot.R
#'
#' @method plot pense
#' @export
plot.pense <- function(x, what = c("cv", "coef.path"), ...) {
    what <- match.arg(what)

    if (what == "cv") {
        .plot_cv_res(x$cv_lambda_grid, x$lambda_opt)
    } else {
        .plot_coef_path(x)
    }
}

#' Plot Method for Cross-Validated Elastic Net Models
#'
#' Plot the cross-validation error or the coefficient path for a fitted
#' elastic net regression model.
#'
#' @param x a fitted, cross-validated EN model from \code{\link{elnet_cv}}.
#' @param what plot either the cross-validated prediction error
#'      (\code{"cv"}; default) or the coefficient paths.
#' @param ... currently ignored.
#'
#' @example examples/elnet_cv-plot.R
#'
#' @method plot cv_elnetfit
#'
#' @export
plot.cv_elnetfit <- function(x, what = c("cv", "coef.path"), ...) {
    what <- match.arg(what)

    if (what == "cv") {
        .plot_cv_res(x$cvres, x$lambda_opt)
    } else {
        .plot_coef_path(x)
    }
}

#' Plot Method for Fitted Elastic Net Models
#'
#' Plot the coefficient path for a fitted elastic net regression model.
#'
#' @param x a fitted EN model from \code{\link{elnet}}.
#' @param ... currently ignored.
#'
#' @example examples/elnet-plot.R
#'
#' @method plot elnetfit
#' @export
plot.elnetfit <- function(x, ...) {
    .plot_coef_path(x)
}


#' @importFrom graphics plot
.plot_cv_res <- function(cv_grid, lambda_opt) {
    plot(
        cvavg ~ lambda,
        cv_grid,
        log = "xy",
        col = 1L + (cv_grid$lambda == lambda_opt),
        pch = 20L,
        xlab = expression(lambda),
        ylab = "Robust RMSPE",
        main = "CV Prediction Error"
    )
}

#' @importFrom graphics plot lines text abline
.plot_coef_path <- function(x) {
    var_names <- if (is.null(rownames(x$coefficients))) {
        sprintf("%d", seq_len(nrow(x$coefficients) - 1L))
    } else {
        rownames(x$coefficients)[-1L]
    }

    active_vars <- data.frame(
        var = var_names[x$coefficients@i],
        lambda = rep.int(x$lambda, times = diff(x$coefficients@p) - 1L),
        value = x$coefficients@x[x$coefficients@i > 0L]
    )

    var_names <- levels(active_vars$var)
    var_colors <- (seq_along(var_names) - 1L) %% 10L + 1L
    names(var_colors) <- var_names

    plot(
        1, 0,
        xlim = range(x$lambda),
        ylim = range(0, active_vars$value),
        type = "n",
        log = "x",
        xlab = expression(lambda),
        ylab = "Coefficient"
    )
    for (i in seq_along(var_names)) {
        sel_var <- levels(active_vars$var)[[i]]
        with(active_vars, {
            act_values <- value[var == sel_var]
            act_lambda <- lambda[var == sel_var]
            matching_lambdas <- match(act_lambda, x$lambda)
            lambda_bound <- range(matching_lambdas, na.rm = TRUE)
            pad_0s <- 0L
            if (lambda_bound[1L] > 1L) {
                act_lambda <- c(
                    act_lambda,
                    x$lambda[lambda_bound[[1L]] - 1L],
                    min(lambda)
                )
                pad_0s <- 2L
            }
            if (lambda_bound[2L] < length(x$lambda)) {
                act_lambda <- c(
                    act_lambda,
                    x$lambda[lambda_bound[[2L]] + 1L],
                    max(lambda)
                )
                pad_0s <- pad_0s + 2L
            }
            lines(
                x = act_lambda,
                y = c(act_values, numeric(pad_0s)),
                col = var_colors[[sel_var]]
            )
            text(
                x = min(act_lambda),
                y = act_values[which.min(act_lambda)],
                labels = sel_var,
                adj = 1,
                pos = 2,
                cex = 0.9,
                col = "#4f4f4f"
            )
        })
    }

    if (!is.null(x$lambda_opt)) {
        abline(v = x$lambda_opt, lty = 2, col = "#4f4f4f")
    }
}

