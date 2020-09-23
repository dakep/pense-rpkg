library(pense)
library(testthat)

test_that("Using a parallel cluster", {
  library(parallel)
  cl <- makePSOCKcluster(1)
  on.exit(stopCluster(cl), add = TRUE, after = FALSE)

  # User messes with cluster...
  clusterEvalQ(cl, {
    x <- "x"
    y <- "y"
    std_data <- "standardized data"
  })

  # Generate dummy data with n=50 observations and p=25 possible predictors
  # (of which only the first 3 are truly relevant).
  n <- 50
  p <- 20
  set.seed(123)
  x <- matrix(rt(n * p, df = 5), ncol = p)
  y <- x[, 1] + 0.5 * x[, 2] + 2 * x[, 3] + rt(n, df = 2)

  set.seed(123) # Setting the seed is suggested for reproducibility of the CV results.
  fit_with_cluster <- adapense_cv(x, y, nlambda = 25, alpha = 0.9, cv_k = 3, cv_repl = 2, cl = cl)
  expect_is(fit_with_cluster, 'pense_cvfit')
})
