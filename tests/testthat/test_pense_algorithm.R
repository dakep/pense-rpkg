library(pense)
library(testthat)

test_that("PENSE Algorithm (1 thread)", {
  skip_if_not_installed('jsonlite')
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')

  requireNamespace('jsonlite')

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

  snapshot_file <- test_path('snap/pense_algo_st.json')
  if (!file.exists(snapshot_file)) {
    jsonlite::write_json(pr$estimates, path = snapshot_file, auto_unbox = TRUE,
                         digits = 14, pretty = TRUE)
    skip('Snapshot file did not exist and was created.')
  }

  ref_estimates <- jsonlite::read_json(snapshot_file, simplifyVector = FALSE)

  for (i in seq_along(ref_estimates)) {
    for (name in names(ref_estimates[[i]])) {
      expect_success(expect_equal(drop(pr$estimates[[!!i]][[!!name]]),
                                  unlist(ref_estimates[[!!i]][[!!name]]),
                                  tolerance = 1e-6))
    }
  }
})

test_that("PENSE Algorithm (2 threads)", {
  skip_if_not_installed('jsonlite')
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')
  skip_if_not(pense:::.k_multithreading_support,
              message = "System does not support multithreading.")

  requireNamespace('jsonlite')

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

  snapshot_file <- test_path('snap/pense_algo_mt.json')
  if (!file.exists(snapshot_file)) {
    jsonlite::write_json(pr$estimates, path = snapshot_file, auto_unbox = TRUE,
                         digits = 14, pretty = TRUE)
    skip('Snapshot file did not exist and was created.')
  }

  ref_estimates <- jsonlite::read_json(snapshot_file, simplifyVector = FALSE)

  for (i in seq_along(ref_estimates)) {
    for (name in names(ref_estimates[[i]])) {
      expect_success(expect_equal(drop(pr$estimates[[!!i]][[!!name]]),
                                  unlist(ref_estimates[[!!i]][[!!name]]),
                                  tolerance = 1e-3))
    }
  }
})

