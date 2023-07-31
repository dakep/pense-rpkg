library(pense)
library(testthat)

test_that("PENSE Algorithm (1 thread)", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 80L
  p <- 20L
  sparse <- FALSE
  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)
  pr <- pense(x, y, alpha = 0.8, nlambda = nlambda, nlambda_enpy = 5,
              ncores = 1L, bdp = 0.25,
              sparse = sparse, eps = 1e-8,
              enpy_opts = enpy_options(retain_max = 5,
                                       en_algorithm_opts = en_lars_options()),
              algorithm_opts = mm_algorithm_options(
                en_algorithm_opts = en_lars_options()))

  expect_equal(pr$bdp, 0.2469)
  expect_equal(pr$alpha, 0.8)
  expect_length(pr$estimates, nlambda)

  compare_estimates(pr$estimates, 'snap/pense_algo_st.json', tol = 1e-6)
})

test_that("PENSE Algorithm (2 threads)", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')
  skip_if_not(pense:::.k_multithreading_support,
              message = "System does not support multithreading.")

  n <- 80L
  p <- 40L
  sparse <- FALSE
  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)
  pr <- pense(x, y, alpha = 0.8,
              nlambda = nlambda, nlambda_enpy = 5,
              ncores = 2L, bdp = 0.25,
              sparse = sparse, eps = 1e-8,
              enpy_opts = enpy_options(retain_max = 5,
                                       en_algorithm_opts = en_lars_options()),
              algorithm_opts = mm_algorithm_options(
                en_algorithm_opts = en_lars_options()))

  expect_equal(pr$bdp, 0.2469)
  expect_equal(pr$alpha, 0.8)
  expect_length(pr$estimates, nlambda)

  compare_estimates(pr$estimates, 'snap/pense_algo_mt.json', tol = 1e-3)
})

