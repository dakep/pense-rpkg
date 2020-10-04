## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", error = TRUE, warning = TRUE)
load("computing_adapense_fits.RData")

## -----------------------------------------------------------------------------
library(pense)

## -----------------------------------------------------------------------------
library(parallel)
# If you don't know how many CPU cores are available, first run `detectCores(logical = FALSE)`
cluster <- makeCluster(3)

## -----------------------------------------------------------------------------
set.seed(1234)
x <- matrix(rweibull(50 * 40, 2, 3), ncol = 40)
y <- 1 * x[, 1] - 1.1 * x[, 2] + 1.2 * x[, 3] + rt(nrow(x), df = 2)

## -----------------------------------------------------------------------------
y[1:3] <- 5 * apply(x[1:3, ], 1, max)
x[3:6, 4:6] <- 1.5 * max(x) + abs(rcauchy(4 * 3))

## ----fit_075, eval=FALSE------------------------------------------------------
#  set.seed(1234)
#  fit_075 <- adapense_cv(x, y, alpha = 0.75, cv_k = 5, cv_repl = 3, cl = cluster)

## ---- include=FALSE-----------------------------------------------------------
fit_075 <- fit_075_3repl

## ---- fig.width=7, fig.height=5, fig.cap="Estimated prediction accuracy using 3 replications of 5-fold CV."----
plot(fit_075)

## ---- include=FALSE-----------------------------------------------------------
fit_075 <- fit_075_10repl

## ---- eval=FALSE--------------------------------------------------------------
#  set.seed(1234)
#  fit_075 <- adapense_cv(x, y, alpha = 0.75, cv_k = 5, cv_repl = 10, cl = cluster)
#  plot(fit_075)

## ---- fig.width=7, fig.height=5, echo=FALSE, fig.cap="Estimated prediction accuracy using 10 replications of 5-fold CV."----
plot(fit_075)

## -----------------------------------------------------------------------------
summary(fit_075)

## -----------------------------------------------------------------------------
summary(fit_075, lambda = "se")

## ----fit_all, eval=FALSE------------------------------------------------------
#  set.seed(1234)
#  fit_075_2 <- adapense_cv(x, y, alpha = 0.75, exponent = 2, cv_k = 5, cv_repl = 10, cl = cluster)
#  
#  set.seed(1234)
#  fit_100_1 <- adapense_cv(x, y, alpha = 1, exponent = 1, cv_k = 5, cv_repl = 10, cl = cluster)
#  
#  set.seed(1234)
#  fit_100_2 <- adapense_cv(x, y, alpha = 1, exponent = 2, cv_k = 5, cv_repl = 10, cl = cluster)

## -----------------------------------------------------------------------------
prediction_performance(fit_075, fit_075_2, fit_100_1, fit_100_2)

## -----------------------------------------------------------------------------
prediction_performance(fit_075, fit_075_2, fit_100_1, fit_100_2, lambda = 'se')

## -----------------------------------------------------------------------------
summary(fit_075_2, lambda = 'se')

## ---- eval=FALSE--------------------------------------------------------------
#  mae <- function (prediction_errors) {
#    mean(abs(prediction_errors))
#  }
#  
#  set.seed(1234)
#  fit_075_mae <- adapense_cv(x, y, alpha = 0.75, cv_k = 5, cv_repl = 5, cl = cluster, cv_metric = mae)

