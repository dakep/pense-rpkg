library(pense)
library(testthat)

test_that("RIS-CV", {
  skip_if_not(nzchar(Sys.getenv('PENSE_TEST_FULL')),
              message = 'Environment variable `PENSE_TEST_FULL` not defined.')
  n <- 50L
  p <- 10L
  nlambda <- 25L
  alphas <- c(0.1, 0.8)

  set.seed(123)
  x <- matrix(rcauchy(n * p), ncol = p)
  y <- 2 + rowSums(x[, 1:5]) / 5 + rnorm(n, sd = 4)

  max_solutions <- 10

  pr <- pense_cv(x, y,
                 cv_k = 3, cv_repl = 2,
                 cv_type = 'ris',
                 alpha = alphas,
                 nlambda = nlambda,
                 nlambda_enpy = 5,
                 max_solutions = max_solutions,
                 ncores = 2L,
                 bdp = 0.25,
                 sparse = FALSE, eps = 1e-8,
                 enpy_opts = enpy_options(retain_max = 5, en_algorithm_opts = en_lars_options()),
                 algorithm_opts = mm_algorithm_options(en_algorithm_opts = en_lars_options()))

  expect_length(pr$estimates, length(alphas) * nlambda)
  expect_contains(colnames(pr$cvres), c('cvavg', 'cvse', 'lambda', 'alpha'))
  expect_length(coef(pr), p + 1L)
})
