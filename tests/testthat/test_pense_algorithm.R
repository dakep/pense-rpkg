library(pense)
library(testthat)

source(test_path('helper.R'))

test_that("PENSE Algorithm (\"best\" solution) (1 thread)", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 80L
  p <- 20L
  sparse <- FALSE
  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)
  pr <- pense(x, y,
              alpha = 0.8,
              nlambda = nlambda,
              nlambda_enpy = 5,
              ncores = 1L,
              bdp = 0.25,
              sparse = sparse,
              eps = 1e-8,
              explore_solutions = 0,
              max_solutions = 1,
              enpy_opts = enpy_options(retain_max = 5,
                                       en_algorithm_opts = en_lars_options()),
              algorithm_opts = mm_algorithm_options(
                en_algorithm_opts = en_lars_options()))

  expect_equal(pr$bdp, 0.2469)
  expect_equal(pr$alpha, 0.8)
  expect_length(pr$estimates, nlambda)

  compare_estimates(pr$estimates, 'snap/pense_algo_st.json', tol = 1e-6)
})

test_that("PENSE Algorithm (multiple solutions)", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 80L
  p <- 20L
  sparse <- FALSE
  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)
  expect_no_condition(
    pr <- pense(x, y,
                alpha = 0.8,
                nlambda = nlambda,
                nlambda_enpy = 5,
                ncores = 1L,
                bdp = 0.25,
                sparse = sparse,
                eps = 1e-8,
                explore_solutions = 2,
                max_solutions = 10,
                enpy_opts = enpy_options(retain_max = 5,
                                         en_algorithm_opts = en_lars_options()),
                algorithm_opts = mm_algorithm_options(
                  en_algorithm_opts = en_lars_options()))
  )

  expect_equal(pr$bdp, 0.2469)
  expect_equal(pr$alpha, 0.8)
  expect_length(pr$estimates, nlambda)
  expect_all_true(lengths(pr$estimates) <= 10)
})

test_that("PENSE Algorithm (\"best\" solution) (2 threads)", {
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

  expect_no_condition(
    pr <- pense(x, y,
                alpha = 0.8,
                nlambda = nlambda,
                nlambda_enpy = 5,
                ncores = 2L,
                bdp = 0.25,
                sparse = sparse,
                eps = 1e-8,
                max_solutions = 1,
                explore_solutions = 0,
                enpy_opts = enpy_options(retain_max = 5,
                                         en_algorithm_opts = en_lars_options()),
                algorithm_opts = mm_algorithm_options(
                  en_algorithm_opts = en_lars_options()))
  )

  expect_equal(pr$bdp, 0.2469)
  expect_equal(pr$alpha, 0.8)
  expect_length(pr$estimates, nlambda)
})


test_that("PENSE Algorithm (optimal rho)", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 80L
  p <- 20L
  sparse <- FALSE
  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)
  expect_no_condition(
    pr <- pense(x, y,
                alpha = 0.8,
                nlambda = nlambda,
                nlambda_enpy = 5,
                ncores = 1L,
                bdp = 0.25,
                sparse = sparse,
                eps = 1e-8,
                explore_solutions = 0,
                max_solutions = 1,
                mscale_opts = mscale_algorithm_options("mopt"),
                enpy_opts = enpy_options(retain_max = 5,
                                         en_algorithm_opts = en_lars_options()),
                algorithm_opts = mm_algorithm_options(
                  en_algorithm_opts = en_lars_options()))
  )

  expect_equal(pr$bdp, 0.2469)
  expect_equal(pr$alpha, 0.8)
  expect_length(pr$estimates, nlambda)

  compare_estimates(pr$estimates, 'snap/pense_algo_st-mopt.json', tol = 1e-4)
})
