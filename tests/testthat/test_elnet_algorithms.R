library(pense)
library(testthat)

## Gradient of the LS-EN objective function
##
## (1/2n) sum(weights * residuals^2) +
##   + lambda * (0.5 * (1 - alpha) * sum(beta^2) + alpha * penalty_loadings * abs(beta))
##
grad_en <- function (est, x, y, weights, intercept = TRUE, penalty_loadings,
                     lambda = est$lambda) {
  if (missing(weights) || is.null(weights)) {
    weights <- rep.int(1, length(y))
  }
  if (missing(penalty_loadings) || is.null(penalty_loadings)) {
    penalty_loadings <- rep.int(1, ncol(x))
  }

  # Loss + L2
  lambda_2 <- lambda * (1 - est$alpha) # 0.5 is remove by square!
  resid <- drop(y - x %*% est$beta - est$intercept)
  grad_beta <- -drop((resid * weights) %*% x) / length(y) +
    as.numeric(lambda_2 * est$beta)
  grad <- if (intercept) {
    c(-mean(resid * weights), grad_beta)
  } else {
    grad_beta
  }
  if (est$alpha < .Machine$double.eps) {
    return(grad)
  }
  # + subgradient
  al <- est$alpha * lambda * penalty_loadings

  nnz_betas <- as.vector(abs(est$beta) > 0)
  if (!is.logical(nnz_betas)) {
    fail("Cannot determine indices of non-zero coefficient")
    return(0)
  }

  nzind <- intercept + which(nnz_betas)
  zind <- intercept + which(abs(est$beta) == 0)
  if (length(nzind) > 0L) {
    grad[nzind] <- grad[nzind] + al[nzind - intercept] * sign(est$beta[nzind - intercept])
  }
  grad[zind] <- ifelse(grad[zind] < -al[zind - intercept],
                       grad[zind] + al[zind - intercept],
                       ifelse(grad[zind] > al[zind - intercept],
                              grad[zind] - al[zind - intercept],
                              0))
  return(grad)
}

check_en_algorithm <- function (en_algorithm_opts, alphas = c(0.1, 0.5, 0.8, 1), num_tol, num_tol_comp = num_tol) {
  set.seed(123)
  x <- matrix(rnorm(100 * 10), ncol = 10)
  y <- 6 + rowSums(x[ , 1:3]) + rnorm(nrow(x))
  obs_wgts <- runif(length(y), 0.5, 1)
  pen_loadings <- runif(ncol(x), 0.1, 1)
  lambda <- c(0.8, 0.4, 0.1)

  max_gradient <- function (alpha, weights, loadings, intercept) {
    weights <- !grepl('w/o', weights, fixed = TRUE)
    loadings <- !grepl('w/o', loadings, fixed = TRUE)
    intercept <- !grepl('w/o', intercept, fixed = TRUE)
    wgts <- if (weights) {
      obs_wgts
    } else {
      NULL
    }
    ldgs <- if (loadings) {
      pen_loadings
    } else {
      NULL
    }
    ests <- if (loadings) {
      if (weights) {
        elnet(x, y, alpha = alpha, lambda = lambda, weights = wgts, penalty_loadings = ldgs, standardize = FALSE,
              en_algorithm_opts = en_algorithm_opts, intercept = intercept, eps = num_tol)
      } else {
        elnet(x, y, alpha = alpha, lambda = lambda, penalty_loadings = ldgs, standardize = FALSE,
              en_algorithm_opts = en_algorithm_opts, intercept = intercept, eps = num_tol)
      }
    } else {
      if (weights) {
        elnet(x, y, alpha = alpha, lambda = lambda, weights = wgts, standardize = FALSE,
              en_algorithm_opts = en_algorithm_opts, intercept = intercept, eps = num_tol)
      } else {
        elnet(x, y, alpha = alpha, lambda = lambda, standardize = FALSE,
              en_algorithm_opts = en_algorithm_opts, intercept = intercept, eps = num_tol)
      }
    }

    max(sapply(ests$estimates, function (est) {
      sum(grad_en(est, x = x, y = y, weights = wgts, penalty_loadings = ldgs, intercept = intercept)^2)
    }))
  }

  invisible(lapply(alphas, function (alpha) {
    lapply(c('w/o weights', 'w/ weights'), function (weights) {
      lapply(c('w/o intercept', 'w/ intercept'), function (intercept) {
        if (alpha > 0) {
          lapply(c('w/o loadings', 'w/ loadings'), function (loadings) {
            expect_lte(max_gradient(alpha = !!alpha, !!weights, !!loadings, !!intercept), !!num_tol_comp)
          })
        } else {
          expect_lte(max_gradient(alpha = !!alpha, !!weights, 'w/o loadings', !!intercept), !!num_tol_comp)
        }
      })
    })
  }))
}

test_that("Elastic Net Algorithm `DAL`", {
  check_en_algorithm(en_dal_options(), num_tol_comp = 1e-7, num_tol = 1e-9)
})

test_that("Elastic Net Algorithm `linearized ADMM`", {
  check_en_algorithm(en_admm_options(), num_tol_comp = 1e-3, num_tol = 1e-12)
})

test_that("Elastic Net Algorithm `augmented LARS`", {
  check_en_algorithm(en_lars_options(), num_tol = 1e-12)
})

test_that("Ridge Algorithm", {
  check_en_algorithm(NULL, alphas = 0, num_tol = 1e-12)
})
