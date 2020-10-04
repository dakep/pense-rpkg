## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", error = TRUE, warning = TRUE)
load("lambda_grids_fits.RData")

## -----------------------------------------------------------------------------
library(pense)

## -----------------------------------------------------------------------------
set.seed(1234)
x <- matrix(rweibull(50 * 40, 2, 3), ncol = 40)
y <- 1 * x[, 1] - 1.1 * x[, 2] + 1.2 * x[, 3] + rt(nrow(x), df = 2)

## -----------------------------------------------------------------------------
y[1:3] <- 5 * apply(x[1:3, ], 1, max)
x[3:6, 4:6] <- 1.5 * max(x) + abs(rcauchy(4 * 3))

## ---- eval=FALSE--------------------------------------------------------------
#  set.seed(1234)
#  fit_grid_narrow <- adapense_cv(x, y, alpha = 0.75, lambda_min_ratio = 1e-1, cv_k = 5, cv_repl = 10)

## ---- eval=FALSE--------------------------------------------------------------
#  set.seed(1234)
#  fit_grid_wide <- adapense_cv(x, y, alpha = 0.75, lambda_min_ratio = 1e-6, cv_k = 5, cv_repl = 10)

## ---- fig.width=6.5, fig.show='hold', echo=FALSE, fig.cap="Prediction performance of models estimated on different grids of the penalization level: (a) narrow grid with `lambda_min_ratio=1e-1`, (b) wide grid with `lambda_min_ratio=1e-6`."----
layout(matrix(1:2, nrow = 1, byrow = TRUE))
prev_par <- par(mai = c(0.8, 0.8, 0.6, 0.2), oma = c(0, 0, 0, 0), lwd = 1.6, cex.main = 0.9, mgp = c(2.5, 1, 0))

plot(fit_grid_narrow)
title(main = '(a) narrow grid\n\n', cex.main = 1)

plot(fit_grid_wide)
title(main = '(b) wide grid\n\n', cex.main = 1)
par(prev_par)

## ---- eval=FALSE--------------------------------------------------------------
#  set.seed(1234)
#  fit_grid_focused <- adapense_cv(x, y, alpha = 0.75, lambda_min_ratio = 1e-2, cv_k = 5, cv_repl = 10)

## ---- echo=FALSE, fig.width=5, fig.height=4, fig.cap="Prediction performance of models estimated over a well-focused grid of penalization levels."----
plot(fit_grid_focused)

## ---- eval=FALSE--------------------------------------------------------------
#  fit_preliminary <- pense_cv(x, y, alpha = 0, cv_k = 5, cv_repl = 10, lambda = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1))
#  exponent <- 1
#  penalty_loadings <- 1 / abs(coef(fit_preliminary)[-1])^exponent
#  fit_adaptive <- pense_cv(x, y, alpha = 0.75, cv_k = 5, cv_repl = 10, lambda = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1, 5))

