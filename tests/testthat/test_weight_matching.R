library(pense)
library(testthat)

test_that("Weight matching", {
  CV_K <- 3
  CV_REPL <- 5
  N <- 50
  N_SOL <- 4L
  N_LAMBDA <- 6L

  set.seed(123)
  r <- array(rnorm(N_SOL * N * N_LAMBDA),
             dim = c(N, N_SOL, N_LAMBDA))

  global_ests <- apply(r, 3, simplify = FALSE, FUN = function (rlambda) {
    apply(rlambda, 2, simplify = FALSE, FUN = function (rsol) {
      list(residuals = rsol, scale = abs(rsol[[1]]))
    })
  })

  set.seed(1234)
  cv_estimates <- lapply(seq_len(CV_K * CV_REPL), function (.) {
    test_ind <- sort(sample.int(N, floor(N / CV_K)))
    train_ind <- seq_len(N)[-test_ind]

    estimates <- apply(r, 3, simplify = FALSE, FUN = function (rlambda) {
      apply(rlambda, 2, simplify = FALSE, FUN = function (rsol) {
        list(residuals = rsol[train_ind],
             test_residuals = rsol[test_ind],
             scale = abs(rsol[[1]]))
      })
    })

    list(test_ind = test_ind,
         train_ind = train_ind,
         estimates = estimates)
  })

  matches <- .Call(pense:::C_match_solutions_by_weight,
                   cv_estimates,
                   global_ests,
                   as.integer(CV_K),
                   as.integer(N),
                   2.5, # cc
                   1L)  # nthreads

  expect_length(matches, N_LAMBDA)

  for (lambda_ind in seq_along(matches)) {
    expect_length(matches[[!!lambda_ind]], N_SOL)
    for (sol_ind in seq_along(matches[[lambda_ind]])) {
      expect_equal(matches[[!!lambda_ind]][[!!sol_ind]]$rankcorr,
                   matrix(1, nrow = CV_K, ncol = CV_REPL))
      expect_length(matches[[!!lambda_ind]][[!!sol_ind]]$wmspe, CV_REPL)
    }
  }

  if (isTRUE(pense:::.k_multithreading_support)) {
    matches_mt <- .Call(pense:::C_match_solutions_by_weight,
                        cv_estimates,
                        global_ests,
                        as.integer(CV_K),
                        as.integer(N),
                        2.5, # cc
                        2L)  # nthreads

    expect_equal(matches_mt, matches)
  }

  # Reverse order of CV solutions
  cv_estimates_rev <- lapply(cv_estimates, function (cve) {
    cve$estimates <- lapply(cve$estimates, rev)
    cve
  })

  matches_rev <- .Call(pense:::C_match_solutions_by_weight,
                       cv_estimates_rev,
                       global_ests,
                       as.integer(CV_K),
                       as.integer(N),
                       2.5, # cc
                       1L)  # nthreads
  expect_equal(matches_mt, matches_rev)
})


test_that("Weight matching (single lambda/solution)", {
  CV_K <- 3
  CV_REPL <- 5
  N <- 50
  N_SOL <- 1L
  N_LAMBDA <- 1L

  set.seed(123)
  r <- array(rnorm(N_SOL * N * N_LAMBDA),
             dim = c(N, N_SOL, N_LAMBDA))

  global_ests <- apply(r, 3, simplify = FALSE, FUN = function (rlambda) {
    apply(rlambda, 2, simplify = FALSE, FUN = function (rsol) {
      list(residuals = rsol, scale = abs(rsol[[1]]))
    })
  })

  set.seed(1234)
  cv_estimates <- lapply(seq_len(CV_K * CV_REPL), function (.) {
    test_ind <- sort(sample.int(N, floor(N / CV_K)))
    train_ind <- seq_len(N)[-test_ind]

    estimates <- apply(r, 3, simplify = FALSE, FUN = function (rlambda) {
      apply(rlambda, 2, simplify = FALSE, FUN = function (rsol) {
        list(residuals = rsol[train_ind],
             test_residuals = rsol[test_ind],
             scale = abs(rsol[[1]]))
      })
    })

    list(test_ind = test_ind,
         train_ind = train_ind,
         estimates = estimates)
  })

  matches <- .Call(pense:::C_match_solutions_by_weight,
                   cv_estimates,
                   global_ests,
                   as.integer(CV_K),
                   as.integer(N),
                   2.5, # cc
                   1L)  # nthreads

  expect_length(matches, N_LAMBDA)

  for (lambda_ind in seq_along(matches)) {
    expect_length(matches[[!!lambda_ind]], N_SOL)

    for (sol_ind in seq_along(matches[[lambda_ind]])) {
      expect_equal(matches[[!!lambda_ind]][[!!sol_ind]]$rankcorr,
                   matrix(1, nrow = CV_K, ncol = CV_REPL))
      expect_length(matches[[!!lambda_ind]][[!!sol_ind]]$wmspe, CV_REPL)
    }
  }

  if (isTRUE(pense:::.k_multithreading_support)) {
    matches_mt <- .Call(pense:::C_match_solutions_by_weight,
                        cv_estimates,
                        global_ests,
                        as.integer(CV_K),
                        as.integer(N),
                        2.5, # cc
                        2L)  # nthreads

    expect_equal(matches_mt, matches)
  }
})
