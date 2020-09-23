library(pense)
library(testthat)

test_that("Response with no variation", {
  # Generate dummy data with n=50 observations and p=25 possible predictors
  # (of which only the first 3 are truly relevant).
  n <- 50
  p <- 5
  set.seed(123)
  x <- matrix(rt(n * p, df = 5), ncol = p)
  y <- c(rep.int(0, n / 2 - 1), rep.int(1, n / 2 + 1)) # has a MAD of 0

  fit_cont <- expect_warning(pense_cv(x, y, alpha = 0.9, nlambda = 10, cv_k = 3), regexp = 'binary')
  expect_is(fit_cont, 'pense_cvfit')
  expect_equal(fit_cont$cv_measure, 'tau_size')
  fit_bin <- expect_warning(pense_cv(x, factor(y), alpha = 0.9, nlambda = 10, cv_k = 3),
                            regexp = 'Cannot compute M-estimate of location')
  expect_is(fit_bin, 'pense_cvfit')
  expect_equal(fit_bin$cv_measure, 'auroc')
})
