## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", error = TRUE, warning = TRUE)

## -----------------------------------------------------------------------------
library(pense)
packageVersion("pense")

## -----------------------------------------------------------------------------
set.seed(1234)
x <- matrix(rt(50 * 10, df = 5), ncol = 10)
y <- 0.5 * x[, 1] - 2 * x[, 2] + 1.5 * x[, 3] + rt(nrow(x), df = 3)

## -----------------------------------------------------------------------------
set.seed(1234)
fitted_with_cv <- pense(x, y, alpha = 0.6, nlambda = 40, warm_reset = 5L, cv_k = 5)

## ---- cache=TRUE--------------------------------------------------------------
set.seed(1234)
fitted_with_cv <- pense_cv(x, y, alpha = 0.6, nlambda = 40, warm_reset = 5L, cv_k = 5)

## ---- cache=TRUE--------------------------------------------------------------
set.seed(1234)
fitted_with_cv <- pense_cv(x, y, alpha = 0.6, nlambda = 40, nlambda_enpy = 5L, cv_k = 5)

## ---- cache=TRUE--------------------------------------------------------------
fitted_no_cv <- pense(x, y, alpha = 0.6, nlambda = 40, nlambda_enpy = 5L)

## -----------------------------------------------------------------------------
str(fitted_no_cv, max.level = 1)
str(fitted_with_cv, max.level = 1)

## -----------------------------------------------------------------------------
coefficients(fitted_with_cv)

## -----------------------------------------------------------------------------
coefficients(fitted_no_cv, lambda = fitted_no_cv$lambda[10])

## -----------------------------------------------------------------------------
coefficients(fitted_with_cv, correction = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  pense(x, y, alpha = 0.6, nlambda = 40, nlambda_enpy = 5L, options = pense_options(delta = 0.33))

## ---- eval=FALSE--------------------------------------------------------------
#  pense(x, y, alpha = 0.6, nlambda = 40, nlambda_enpy = 5L, bdp = 0.33)

## ---- eval=FALSE--------------------------------------------------------------
#  pense(x, y, alpha = 0.6, nlambda = 40, nlambda_enpy = 5L,
#        algorithm_opts = mm_algorithm_options(en_algorithm_opts = en_lars_options()),
#        enpy_opts = enpy_options(en_algorithm_opts = en_admm_options()))

