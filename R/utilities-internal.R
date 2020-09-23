## Get the Constant for Consistency for the M-Scale Using the Bisquare Rho Function
## @param delta desired breakdown point (between 0 and 0.5)
##
## @return consistency constant
#' @importFrom stats pnorm uniroot
.bisquare_consistency_const <- function (delta) {
  ##
  ## Pre-computed values for some delta values
  ##
  eps <- sqrt(.Machine$double.eps)
  if (!isTRUE(delta < 0.5 + eps && delta > -eps)) {
    stop("`delta` is outside valid bounds")
  }

  if (abs(delta - 0.5) < eps) {
    return(1.5476450)
  } else if (abs(delta - 0.25) < eps) {
    return(2.937015)
  } else if (abs(delta - 0.1) < eps) {
    return(5.182361)
  } else if (delta < 0.005) {
    return(50) # ~.1% bdp for bisquare
  }
  integral_interval <- if (delta > 0.1) {
    c(1.5, 5.5)
  } else {
    c(5, 25)
  }

  # For bisquare we have the closed form solution to the expectation
  expectation <- function(cc, delta) {
    pnorm.mcc <- 2 * pnorm(-cc)
    1/cc^6 * exp(-(cc^2/2)) * (
      -cc * (15 - 4 * cc^2 + cc^4) * sqrt(2 / pi) +
        3 * (5 - 3 * cc^2 + cc^4) * exp(cc^2/2) * (1 - pnorm.mcc) +
        cc^6 * exp(cc^2/2) * pnorm.mcc
    ) - delta
  }
  uniroot(expectation, interval = integral_interval, delta)$root
}

## Approximate Value Matching
##
## @param x,table see [base::match] for details.
## @param eps numerical tolerance for matching.
## @return a vector the same length as `x` with integers giving the position in
##         `table` of the first match if there is a match, or `NA_integer_`
##         otherwise.
.approx_match <- function(x, table, eps = min(sqrt(.Machine$double.eps), 0.5 * min(x, table))) {
  .Call(C_approx_match, as.numeric(x), as.numeric(table), as.numeric(eps[[1L]]))
}

## Extract the given metric from all matching nodes (by name).
extract_metric <- function (metrics, attr, node) {
  matches <- c()
  if (!is.null(metrics[[attr]]) && isTRUE(metrics$name == node)) {
    matches <- c(matches, metrics[[attr]])
  }
  if (!is.null(metrics$sub_metrics)) {
    matches <- c(matches, unlist(lapply(metrics$sub_metrics, extract_metric,
                                        attr, node),
                                 use.names = FALSE, recursive = FALSE))
  }
  return (matches)
}

.recurisve_metrics_class <- function (metrics) {
  class(metrics) <- 'nsoptim_metrics'
  if (!is.null(metrics$sub_metrics)) {
    metrics$sub_metrics <- lapply(metrics$sub_metrics, .recurisve_metrics_class)
  }
  return(metrics)
}

.metrics_attrib <- function (estimates, metrics) {
  if (!is.null(metrics) && isTRUE(metrics$name != '')) {
    attr(estimates, 'metrics') <- .recurisve_metrics_class(metrics)
  }
  return(estimates)
}

## Run replicated K-fold CV with random splits
##
## @param std_data standardized full data set (standardized by `.standardize_data`)
## @param cv_k number of folds per CV split
## @param cv_repl number of CV replications.
## @param cv_est_fun function taking the standardized training set and the indices of the left-out observations and
##                   returns a list of estimates. The function always needs to return the same number of estimates!
## @param metric function taking a vector of prediction errors and returning the scale of the prediction error.
#' @importFrom Matrix drop
#' @importFrom rlang abort
.run_replicated_cv <- function (std_data, cv_k, cv_repl, cv_est_fun, metric, par_cluster = NULL) {
  est_fun <- match.fun(cv_est_fun)
  call_with_errors <- isTRUE(length(formals(metric)) == 1L)

  if (length(std_data$y) / cv_k < 2) {
    abort("`cv_k` must be chosen to have at least 2 observations in each fold.")
  }

  test_segments_list <- lapply(integer(cv_repl), function (repl_id) {
    split(seq_along(std_data$y), sample(rep_len(seq_len(cv_k), length(std_data$y))))
  })
  test_segments <- unlist(test_segments_list, recursive = FALSE, use.names = FALSE)

  cl_handler <- .make_cluster_handler(par_cluster)

  predictions_all <- cl_handler(test_segments, function (test_ind, est_fun) {
    train_x <- std_data$x[-test_ind, , drop = FALSE]
    train_y <- std_data$y[-test_ind]
    test_x <- std_data$x[test_ind, , drop = FALSE]

    train_std <- std_data$cv_standardize(train_x, train_y)
    cv_ests <- est_fun(train_std, test_ind)

    matrix(unlist(lapply(cv_ests, function (est) {
      unstd_est <- train_std$unstandardize_coef(est)
      drop(test_x %*% unstd_est$beta) - unstd_est$intercept
    }), use.names = FALSE, recursive = FALSE), ncol = length(cv_ests))
  }, est_fun = est_fun)

  predictions_all <- split(predictions_all, rep(seq_len(cv_repl), each = cv_k))
  prediction_metrics <- mapply(predictions_all, test_segments_list, FUN = function (predictions, test_inds) {
    obs_order <- sort.list(unlist(test_inds, recursive = FALSE, use.names = FALSE))
    ordered_predictions <- do.call(rbind, predictions)[obs_order, ]

    if (call_with_errors) {
      apply(ordered_predictions - std_data$y, 2, metric)
    } else {
      apply(ordered_predictions, 2, metric, std_data$y)
    }
  })
  matrix(unlist(prediction_metrics, recursive = FALSE, use.names = FALSE), ncol = cv_repl)
}

## Standardize data
##
## @param x predictor matrix. Can also be a list with components `x` and `y`, in which case `y` is ignored.
## @param y response vector.
## @param intercept is an intercept included (i.e., should `y` be centered?)
## @param standardize standardize or not.
## @param robust use robust standardization.
## @param location_rho,location_cc rho function and cutoff supplied to `mlocscale()` and `mloc()`
## @param ... passed on to `mlocscale()`.
## @return a list with the following entries:
#' @importFrom Matrix drop
#' @importFrom methods is
#' @importFrom rlang abort
#' @importFrom stats sd
.standardize_data <- function (x, y, intercept, standardize, robust, sparse, mscale_opts, location_rho = 'bisquare',
                               location_cc = 4.5, target_scale_x = NULL, ...) {
  if (is.list(x) && !is.null(x$x) && !is.null(x$y)) {
    y <- x$y
    x <- x$x
  }

  ret_list <- list(scale_x = rep.int(1, ncol(x)), mux = numeric(ncol(x)), muy = 0, x = x, y = y)

  ## Center data for numerical convenience
  if (isTRUE(intercept)) {
    if (!isTRUE(robust)) {
      ret_list$mux <- colMeans(x)
      ret_list$muy <- mean(y)
    } else {
      ret_list$mux <- apply(x, 2, function (xj) {
        mloc(xj, rho = location_rho, cc = location_cc, opts = mscale_opts)
      })
      ret_list$muy <- mloc(y, rho = location_rho, cc = location_cc, opts = mscale_opts)
      if (!is.finite(ret_list$muy)) {
        # In case the response has more than 50% equal values.
        ret_list$muy <- 0
      }
    }

    ret_list$x <- sweep(x, 2L, ret_list$mux, FUN = `-`, check.margin = FALSE)
    ret_list$y <- y - ret_list$muy
  }

  ## Scale predictors
  if (isTRUE(standardize) || isTRUE(standardize == 'cv_only')) {
    ret_list$scale_x <- if (!isTRUE(robust)) {
      apply(ret_list$x, 2, sd)
    } else {
      locscale <- apply(ret_list$x, 2, function (xj) {
        mlocscale(xj, location_rho = location_rho, location_cc = location_cc, opts = mscale_opts, ...)
      })
      if (isTRUE(intercept)) {
        # Re-center the predictors with the updated centers
        ret_list$mux <- ret_list$mux + locscale[1L, ]
        ret_list$x <- sweep(x, 2L, ret_list$mux, FUN = `-`, check.margin = FALSE)
      }
      locscale[2L, ]
    }

    if (!isTRUE(all(ret_list$scale_x > 0))) {
      abort("Standardization failed. One or more variables in `x` have a scale of 0.")
    }

    if (isTRUE(standardize)) {
      ret_list$x <- if (!is.null(target_scale_x)) {
        sweep(ret_list$x, 2L, target_scale_x / ret_list$scale_x, FUN = `*`, check.margin = FALSE)
      } else {
        sweep(ret_list$x, 2L, ret_list$scale_x, FUN = `/`, check.margin = FALSE)
      }
    }
  }

  # Set the target scale to 1, so that standardizing and unstandardizing works.
  if (is.null(target_scale_x)) {
    target_scale_x <- 1
  }

  ret_list$cv_standardize <- function(x, y) {
    if (is.list(x) && !is.null(x$x) && !is.null(x$y)) {
      y <- x$y
      x <- x$x
    }

    if (isTRUE(standardize == 'cv_only')) {
      # In case of "CV only" standardization, match the original scaling
      .standardize_data(x, y, intercept = intercept, standardize = TRUE, robust = robust, sparse = sparse,
                        location_rho = location_rho, location_cc = location_cc, mscale_opts = mscale_opts,
                        target_scale_x = ret_list$scale_x, ... = ...)
    } else {
      .standardize_data(x, y, intercept = intercept, standardize = standardize, robust = robust, sparse = sparse,
                        location_rho = location_rho, location_cc = location_cc, mscale_opts = mscale_opts, ... = ...)
    }
  }

  ret_list$standardize_coefs <- function(coef_obj) {
    if (is.null(coef_obj)) {
      return(NULL)
    }
    if (is.null(coef_obj$intercept)) {
      coef_obj$intercept <- 0
    }
    if (isTRUE(intercept)) {
      # Adjust intercept
      coef_obj$intercept <- coef_obj$intercept - ret_list$muy + sum(ret_list$mux * coef_obj$beta)
    }
    if (isTRUE(standardize)) {
      coef_obj$beta <- coef_obj$beta * (ret_list$scale_x / target_scale_x)
    }
    return(coef_obj)
  }

  ret_list$unstandardize_coefs <- if (isTRUE(sparse)) {
    function(coef_obj) {
      if (is.null(coef_obj)) {
        return(coef_obj)
      }
      if (is.null(coef_obj$intercept)) {
        coef_obj$intercept <- 0
      }
      if (isTRUE(standardize)) {
        coef_obj$beta@x <- coef_obj$beta@x * target_scale_x / ret_list$scale_x[coef_obj$beta@i]
      }
      if (isTRUE(intercept)) {
        # Recreate intercept
        coef_obj$intercept <- coef_obj$intercept + ret_list$muy - sum(ret_list$mux[coef_obj$beta@i] * coef_obj$beta@x)
      }
      return(coef_obj)
    }
  } else {
    function(coef_obj) {
      if (is.null(coef_obj)) {
        return(coef_obj)
      }
      if (is.null(coef_obj$intercept)) {
        coef_obj$intercept <- 0
      }
      if (isTRUE(standardize)) {
        coef_obj$beta <- coef_obj$beta * target_scale_x / ret_list$scale_x
      }
      if (isTRUE(intercept)) {
        # Recreate intercept
        coef_obj$intercept <- coef_obj$intercept + ret_list$muy - sum(ret_list$mux * coef_obj$beta)
      }
      return(coef_obj)
    }
  }

  return(ret_list)
}

#' @importFrom stats median
.cv_mape <- function (r) {
  median(abs(r))
}

.cv_rmspe <- function (r) {
  sqrt(mean(r^2))
}

## Area under the ROC for "negatives" having value 0 and "positives" having value 1.
.cv_auroc <- function (pred, truth) {
  n_neg <- sum(truth <= 0)
  n_pos <- sum(truth > 0)
  mww <- sum(rank(pred)[truth <= 0]) - n_neg * (n_neg + 1) / 2
  return(mww / (n_neg * n_pos))
}

.cv_se_selection <- function (cvm, cvsd, se_fact) {
  type <- rep.int(factor('none', levels = c('none', 'min', 'se_fact')), length(cvm))
  best <- which.min(cvm)
  candidates <- which(cvm <= cvm[[best]] + se_fact * cvsd[[best]])
  candidates <- candidates[candidates <= best]  # only consider sparser solutions

  # "ignore" solutions after which the prediction performance comes back down
  best_1se <- if (any(diff(candidates) > 1)) {
    min(candidates[seq_len(max(which(diff(candidates) > 1)))])
  } else {
    min(candidates)
  }
  type[[best]] <- 'min'
  type[[best_1se]] <- 'se_fact'
  return(type)
}

## Validate the response type and return a list with elements `values` and `binary`.
#' @importFrom rlang warn abort
.validate_response <- function (y) {
  nlevels <- length(unique(y))
  if (is.factor(y)) {
    if (nlevels == 2L) {
      return(list(values = as.numeric(y) - 1, binary = TRUE))
    } else if (nlevels > 2L) {
      warn("`y` with more than 2 factor levels is implicitly treated as numeric.")
      return(list(values = as.numeric(y), binary = FALSE))
    }
  } else {
    y <- .as(y, 'numeric')
    if (nlevels == 2L) {
      warn(paste("`y` is interpreted as continuous response but has only 2 distinct values.",
                 "If binary classification is desired, coerce `y` to factor with `as.factor()`."))
    }
    if (nlevels > 1L) {
      return(list(values = y, binary = FALSE))
    }
  }
  abort("`y` must have at least 2 distinct values.")
}

## A wrapper around `methods::as` which raises an error if the conversion results in NA.
## @param ... passed on to [methods::as].
##
#' @importFrom methods as
#' @importFrom rlang abort
.as <- function (object, class, ...) {
  object_var <- deparse(substitute(object))
  tryCatch(methods::as(object, class, ...), warning = function (w) {
    abort(sprintf('`%s` is not of type `%s`', object_var, class))
  })
}

## Create a function which restores the original length of the coefficient vector
.restore_coef_length_fun <- function (positions, length) {
  if (length(positions) > 0L) {
    function (coef) {
      if (is(coef$beta, 'dsparseVector')) {
        coef$beta <- sparseVector(coef$beta@x, positions[coef$beta@i], length)
      } else {
        beta <- numeric(length)
        beta[positions] <- coef$beta
        coef$beta <- beta
      }
      return(coef)
    }
  } else {
    function (coef) {
      if (is(coef$beta, 'dsparseVector')) {
        coef$beta <- sparseVector(numeric(0L), integer(0L), length)
      } else {
        coef$beta <- numeric(length)
      }
      return(coef)
    }
  }
}

.prepare_penalty_loadings <- function (penalty_loadings, x, alpha, sparse, stop_all_infinite = FALSE) {
  orig_p <- ncol(x)
  restore_coef_length <- function (coef) coef

  if(alpha < .Machine$double.eps) {
    abort("Non-empty `penalty_loadings` only supported for `alpha` > 0.")
  } else if (length(penalty_loadings) != orig_p) {
    abort("`penalty_loadings` has different number of elements than `x` columns.")
  }
  penalty_loadings <- .as(penalty_loadings, 'numeric')

  if (any(penalty_loadings < 0)) {
    abort("`penalty_loadings` must be positive.")
  }

  # Determine finite penalty loadings
  good_pl <- which(is.finite(penalty_loadings))
  if (length(good_pl) < orig_p) {
    # Some penalty loadings are infinite! Remove the corresponding predictors from `x`.
    x <- x[ , good_pl, drop = FALSE]
    penalty_loadings <- penalty_loadings[good_pl]
    restore_coef_length <- if (length(good_pl) > 0L) {
      if (isTRUE(sparse)) {
        function (coef) {
          coef$beta <- sparseVector(coef$beta@x, good_pl[coef$beta@i], orig_p)
          return(coef)
        }
      } else {
        function (coef) {
          beta <- numeric(orig_p)
          beta[good_pl] <- coef$beta
          coef$beta <- beta
          return(coef)
        }
      }
    } else {
      if (stop_all_infinite) {
        abort("At least one value in `penalty_loadings` must be finite.")
      }
      if (isTRUE(sparse)) {
        function (coef) {
          coef$beta <- sparseVector(numeric(0L), integer(0L), orig_p)
          return(coef)
        }
      } else {
        function (coef) {
          coef$beta <- numeric(orig_p)
          return(coef)
        }
      }
    }
  }

  return(list(loadings = penalty_loadings, trimmed_x = x, restore_fun = restore_coef_length))
}


#' @importFrom methods is
.sparsify_other_starts <- function (other_starts, sparse) {
  lapply(other_starts, function (est) {
    if (isTRUE(sparse) && !is(est$beta, 'dsparseVector')) {
      est$beta <- sparseVector(as.numeric(est$beta), seq_along(est$beta), length(est$beta))
    } else if (!isTRUE(sparse) && !is.numeric(est$beta)) {
      est$beta <- .as(est$beta, 'numeric')
    }
    class(est) <- NULL
    return(est)
  })
}


#' @importFrom parallel clusterEvalQ clusterCall clusterApplyLB
#' @importFrom rlang abort
.make_cluster_handler <- function (par_cluster) {
  if (is.null(par_cluster)) {
    return(function (X, FUN, ..., x, fun) {
      lapply(X, FUN, ...)
    })
  } else {
    tryCatch({
      clusterEvalQ(par_cluster, {
        library(pense)
      })
    }, error = function (e) {
      abort(paste("`parallel` cluster cannot be used:", e))
    })

    return(function (X, FUN, ..., x, fun) {
      clusterApplyLB(par_cluster, x = X, fun = function (X, FUN, ...) {
        FUN(X, ...)
      }, FUN = FUN, ... = ...)
    })
  }
}
