library(pense)
library(testthat)

test_that("predict(pense) with multiple alpha & solutions", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 50L
  p <- 10L

  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)
  y[1:4] <- 10 * rowSums(x[1:4, 6:7]) + rnorm(4, sd = 0.2)
  y[5:8] <- 10 * rowSums(x[5:8, 6:7]) + rnorm(4, sd = 0.2)

  pr <- pense_cv(x, y,
                 cv_k = 2,
                 cv_repl = 2,
                 alpha = c(0.1, 0.8),
                 nlambda = nlambda,
                 ncores = 2L,
                 bdp = 0.25,
                 sparse = FALSE, eps = 1e-8)

  expect_equal(pr$alpha, c(0.1, 0.8))
  expect_type(pr$lambda, 'list')
  expect_type(pr$lambda[[1L]], 'double')
  expect_type(pr$lambda[[2L]], 'double')
  expect_length(pr$lambda[[1L]], nlambda)
  expect_length(pr$lambda[[2L]], nlambda)
  expect_length(pr$estimates, 2 * nlambda)

  expect_invisible(plot(pr))
  expect_invisible(plot(pr, alpha = 0.1))
  expect_invisible(plot(pr, alpha = 0.8))
  expect_error(plot(pr, alpha = 0.3), regexp = "available")

  expect_type(coef(pr, lambda = pr$lambda[[1]][[5]], alpha = 0.1), 'double')
  expect_type(coef(pr, lambda = pr$lambda[[1]][[5]], alpha = 0.8), 'double')
  expect_warning(expect_type(coef(pr, lambda = pr$lambda[[1]][[5]], alpha = c(0.1, 0.8)), 'double'))
  expect_warning(coef(pr, lambda = pr$lambda[[1]][[5]]), regexp = 'Using first value')
})

test_that("predict(pense_cv) with multiple alpha & solutions", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 50L
  p <- 10L

  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)

  max_solutions <- 10

  pr <- pense_cv(x, y,
                 cv_k = 3, cv_repl = 10,
                 fit_all = TRUE,
                 alpha = c(0.1, 0.8),
                 nlambda = nlambda,
                 nlambda_enpy = 5,
                 max_solutions = max_solutions,
                 ncores = 2L,
                 bdp = 0.25,
                 sparse = FALSE, eps = 1e-8,
                 enpy_opts = enpy_options(retain_max = 5, en_algorithm_opts = en_lars_options()),
                 algorithm_opts = mm_algorithm_options(en_algorithm_opts = en_lars_options()))

  expect_equal(pr$alpha, c(0.1, 0.8))
  expect_type(pr$lambda, 'list')
  expect_type(pr$lambda[[1L]], 'double')
  expect_type(pr$lambda[[2L]], 'double')
  expect_length(pr$lambda[[1L]], nlambda)
  expect_length(pr$lambda[[2L]], nlambda)
  expect_true(length(pr$estimates) >= 2 * nlambda)
  expect_true(length(pr$estimates) <= 2 * max_solutions * nlambda)
  expect_length(pr$cvres$lambda, 2 * nlambda)
  expect_length(pr$cvres$alpha, 2 * nlambda)
  expect_length(pr$cvres$cvavg, 2 * nlambda)
  expect_length(pr$cvres$cvse, 2 * nlambda)

  expect_invisible(plot(pr, what = 'cv'))
  expect_invisible(plot(pr, what = 'coef'))
  expect_invisible(plot(pr, what = 'coef', alpha = 0.1))
  expect_invisible(plot(pr, what = 'coef', alpha = 0.8))
  expect_error(plot(pr, what = 'coef', alpha = 0.3), regexp = "available")
  expect_error(plot(pr, what = 'cv', alpha = 0.3), regexp = "available")

  expect_type(coef(pr), 'double')
  expect_type(coef(pr, lambda = 'min'), 'double')
  expect_type(coef(pr, alpha = 0.1), 'double')
  expect_type(coef(pr, alpha = 0.8), 'double')
  expect_error(coef(pr, alpha = 0.3), regexp = "not fit")
})
