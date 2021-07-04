library(pense)
library(testthat)

test_that("pense() with multiple alpha", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 50L
  p <- 10L

  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)

  pr <- pense(x, y,
              alpha = c(0.1, 0.8),
              nlambda = nlambda,
              nlambda_enpy = 5,
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

test_that("pense_cv() with multiple alpha", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 50L
  p <- 10L

  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)

  pr <- pense_cv(x, y,
                 cv_k = 3, cv_repl = 2,
                 fit_all = c('min', '1-se', '2-se'),
                 alpha = c(0.1, 0.8),
                 nlambda = nlambda,
                 nlambda_enpy = 5,
                 ncores = 2L,
                 bdp = 0.25,
                 sparse = FALSE, eps = 1e-8,
                 enpy_opts = enpy_options(retain_max = 5, en_algorithm_opts = en_lars_options()),
                 algorithm_opts = mm_algorithm_options(en_algorithm_opts = en_lars_options()))

  expect_equal(pr$alpha, c(0.1, 0.8))
  expect_type(pr$lambda, 'list')
  expect_type(pr$lambda[[1L]], 'double')
  expect_type(pr$lambda[[2L]], 'double')
  expect_length(pr$lambda[[1L]], 2)
  expect_length(pr$lambda[[2L]], 2)
  expect_length(pr$estimates, 4)
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

test_that("regmest() with multiple alpha", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 50L
  p <- 10L

  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)

  pr <- regmest(x, y,
                scale = 1.11,
                alpha = c(0.1, 0.8),
                nlambda = nlambda,
                ncores = 2L,
                sparse = FALSE, eps = 1e-8,
                algorithm_opts = mm_algorithm_options(en_algorithm_opts = en_lars_options()))

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

test_that("regmest_cv() with multiple alpha", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 50L
  p <- 10L

  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)

  pr <- regmest_cv(x, y,
                   cv_k = 3, cv_repl = 2,
                   fit_all = c('min', '1-se', '2-se'),
                   scale = 1.11,
                   alpha = c(0.1, 0.8),
                   nlambda = nlambda,
                   ncores = 2L,
                   bdp = 0.25,
                   sparse = FALSE, eps = 1e-8,
                   algorithm_opts = mm_algorithm_options(en_algorithm_opts = en_lars_options()))

  expect_equal(pr$alpha, c(0.1, 0.8))
  expect_type(pr$lambda, 'list')
  expect_type(pr$lambda[[1L]], 'double')
  expect_type(pr$lambda[[2L]], 'double')
  expect_length(pr$lambda[[1L]], 2)
  expect_length(pr$lambda[[2L]], 2)
  expect_length(pr$estimates, 4)
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
