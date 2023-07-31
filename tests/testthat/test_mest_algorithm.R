library(pense)
library(testthat)

source(test_path('helper.R'))

test_that("Regularized M-estimation Algorithm (1 thread)", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  n <- 80L
  p <- 20L
  sparse <- FALSE
  nlambda <- 25L

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)
  scale <- 2
  pr <- regmest(x, y, alpha = 0.8,
                nlambda = nlambda,
                scale = scale,
                starting_points = c(
                  starting_point(intercept = 2,
                                 beta = c(rep.int(1, 5), numeric(p - 5L))),
                  starting_point(intercept = 1,
                                 beta = c(rep.int(2, 2), numeric(p - 2L)))),
                ncores = 1L,
                mscale_bdp = 0.25,
                sparse = sparse,
                eps = 1e-8,
                algorithm_opts = mm_algorithm_options(
                  en_algorithm_opts = en_lars_options()))

  expect_equal(pr$alpha, 0.8)
  expect_length(pr$estimates, nlambda)

  compare_estimates(pr$estimates, 'snap/regmest_algo_st.json', tol = 1e-6)
})

test_that("Regularized M-estimation Algorithm (2 threads)", {
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
  scale <- 2
  pr <- regmest(x, y, alpha = 0.8,
                nlambda = nlambda,
                scale = scale,
                starting_points = c(
                  starting_point(intercept = 2,
                                 beta = c(rep.int(1, 5), numeric(p - 5L))),
                  starting_point(intercept = 1,
                                 beta = c(rep.int(2, 2), numeric(p - 2L)))),
                ncores = 2L,
                mscale_bdp = 0.25,
                sparse = sparse,
                eps = 1e-8,
                algorithm_opts = mm_algorithm_options(
                  en_algorithm_opts = en_lars_options()))

  expect_equal(pr$alpha, 0.8)
  expect_length(pr$estimates, nlambda)
  compare_estimates(pr$estimates, 'snap/regmest_algo_mt.json', tol = 1e-6)
})

