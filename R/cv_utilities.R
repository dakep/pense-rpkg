#' Run replicated K-fold CV with random splits
#'
#' @param std_data standardized full data set
#'    (standardized by `.standardize_data`)
#' @param cv_k number of folds per CV split
#' @param cv_repl number of CV replications.
#' @param cv_est_fun function taking the standardized training set and
#'    the indices of the left-out observations and returns a list of estimates.
#'    The function always needs to return the same number of estimates!
#' @param metric function taking a vector of prediction errors and
#'    returning the scale of the prediction error.
#' @param par_cluster parallel cluster to parallelize computations.
#' @param handler_args additional arguments to the handler function.
#' @importFrom Matrix drop
#' @importFrom rlang abort
#' @keywords internal
.run_replicated_cv <- function (std_data, cv_k, cv_repl, cv_est_fun, metric,
                                par_cluster = NULL,
                                handler_args = list()) {
  est_fun <- match.fun(cv_est_fun)
  call_with_errors <- isTRUE(length(formals(metric)) == 1L)

  if (length(std_data$y) / cv_k < 2) {
    abort("`cv_k` must be chosen to have at least 2 observations in each fold.")
  }

  test_segments_list <- lapply(integer(cv_repl), function (repl_id) {
    split(seq_along(std_data$y),
          sample(rep_len(seq_len(cv_k), length(std_data$y))))
  })
  test_segments <- unlist(test_segments_list, recursive = FALSE,
                          use.names = FALSE)

  cl_handler <- .make_cluster_handler(par_cluster)

  predictions_all <- cl_handler(
    test_segments,
    function (test_ind, est_fun, handler_args) {
      train_x <- std_data$x[-test_ind, , drop = FALSE]
      train_y <- std_data$y[-test_ind]
      test_x <- std_data$x[test_ind, , drop = FALSE]

      train_std <- std_data$cv_standardize(train_x, train_y)
      cv_ests <- est_fun(train_std, test_ind, handler_args)

      matrix(unlist(lapply(cv_ests, function (est) {
        unstd_est <- train_std$unstandardize_coef(est)
        drop(test_x %*% unstd_est$beta) - unstd_est$intercept
      }), use.names = FALSE, recursive = FALSE), ncol = length(cv_ests))
    }, est_fun = est_fun, handler_args = handler_args)

  predictions_all <- split(predictions_all, rep(seq_len(cv_repl), each = cv_k))
  prediction_metrics <- mapply(
    predictions_all, test_segments_list,
    FUN = function (predictions, test_inds) {
      obs_order <- sort.list(unlist(test_inds, recursive = FALSE, use.names = FALSE))
      ordered_predictions <- do.call(rbind, predictions)[obs_order, ]
      if (is.null(dim(ordered_predictions))) {
        dim(ordered_predictions) <- c(length(ordered_predictions), 1L)
      }

      if (call_with_errors) {
        apply(ordered_predictions - std_data$y, 2, metric)
      } else {
        apply(ordered_predictions, 2, metric, std_data$y)
      }
    },
    SIMPLIFY = FALSE)
  matrix(unlist(prediction_metrics, recursive = FALSE, use.names = FALSE), ncol = cv_repl)
}


#' Run replicated K-fold CV with random splits, matching the global estimates
#' to the CV estimates by Kendall's tau-b computed on the robustness weights.
#'
#' @param std_data standardized full data set
#'    (standardized by `.standardize_data`)
#' @param cv_k number of folds per CV split
#' @param cv_repl number of CV replications.
#' @param cv_est_fun function taking the standardized training set and
#'    the indices of the left-out observations and returns a list of estimates.
#'    The function always needs to return the same number of estimates!
#' @param global_ests estimates computed on all observations.
#' @param min_similarity minimum (average) similarity for CV solutions to be considered
#'   (between 0 and 1).
#'   If no CV solution satisfies this lower bound, the best CV solution will be used regardless
#'   of similarity.
#' @param rho_cc consistency constant for Tukey's bisquare rho function.
#' @param par_cluster parallel cluster to parallelize computations.
#' @param handler_args additional arguments to the handler function.
#' @importFrom Matrix drop
#' @importFrom rlang abort
#' @importFrom stats cor
#' @keywords internal
.run_replicated_cv_ris <- function (std_data, cv_k, cv_repl, cv_est_fun,
                                    global_ests,
                                    min_similarity = 0,
                                    par_cluster = NULL,
                                    rho_cc,
                                    handler_args = list()) {
  est_fun <- match.fun(cv_est_fun)
  if (!isTRUE(min_similarity >= 0) || !isTRUE(min_similarity <= 1)) {
    abort("`min_similarity` must be a scalar between 0 and 1.")
  }

  if (length(std_data$y) / cv_k < 2) {
    abort("`cv_k` must be chosen to have at least 2 observations in each fold.")
  }

  test_segments_list <- lapply(integer(cv_repl), \(repl_id) {
    split(seq_along(std_data$y),
          sample(rep_len(seq_len(cv_k), length(std_data$y))))
  })
  test_segments <- unlist(test_segments_list, recursive = FALSE,
                          use.names = FALSE)

  cl_handler <- .make_cluster_handler(par_cluster)

  dispatcher <- function (test_ind, est_fun, std_data, handler_args) {
    train_x <- std_data$x[-test_ind, , drop = FALSE]
    train_y <- std_data$y[-test_ind]
    test_x <- std_data$x[test_ind, , drop = FALSE]

    train_std <- std_data$cv_standardize(train_x, train_y)

    ## Starting points need to be scaled too
    if (!is.null(handler_args$args$optional_args$individual_starts)) {
      handler_args$args$optional_args$individual_starts <- lapply(
        handler_args$args$optional_args$individual_starts,
        \(starts) { lapply(starts, train_std$standardize_coefs) })
    }
    if (!is.null(handler_args$args$optional_args$shared_starts)) {
      handler_args$args$optional_args$shared_starts <- lapply(
        handler_args$args$optional_args$shared_starts,
        \(starts) { lapply(starts, train_std$standardize_coefs) })
    }

    cv_ests <- est_fun(train_std, test_ind, handler_args)

    train_ind <- seq_along(std_data$y)[-test_ind]

    cv_ests <- lapply(cv_ests, \(ests_lambda) {
      lapply(ests_lambda, \(est) {
        unstd_est <- train_std$unstandardize_coef(est)
        unstd_est$test_residuals <- drop(std_data$y[test_ind] -
                                           test_x %*% unstd_est$beta - unstd_est$intercept)

        # Remove coefficients; they're not needed anymore
        unstd_est$beta <- NULL
        unstd_est$intercept <- NULL
        unstd_est$std_intercept <- NULL
        unstd_est$std_beta <- NULL

        unstd_est
      })
    })

    list(estimates = cv_ests,
         test_ind = test_ind,
         train_ind = train_ind)
  }

  # `estimates_all` will be a flat list of length `cv_repl * cv_k`, each with `nlambda` elements:
  # [element 1 - cv_repl * cv_k]:
  #   [lambda 1 - nlambda]
  #     [solution 1 - Q]

  cv_ests <- cl_handler(test_segments,
                        dispatcher,
                        est_fun = est_fun,
                        std_data = std_data,
                        handler_args = handler_args)

  # Match solutions based on robustness weights
  glbl_wgts <- .Call(C_robustness_weights,
                     global_ests,
                     length(std_data$y),
                     rho_cc)

  matches <- lapply(global_ests, \(x) {
    lapply(x, \(...) {
      list(similarity  = matrix(NA_real_, nrow = cv_k, ncol = cv_repl),
           cv_avg_wgt  = matrix(0,        nrow = length(std_data$y), ncol = cv_repl),
           predictions = matrix(NA_real_, nrow = length(std_data$y), ncol = cv_repl))
    })
  })

  for (cv_ind in seq_along(cv_ests)) {
    cvr <- 1L + (cv_ind - 1L) %/% cv_k
    x <- cv_ests[[cv_ind]]
    cv_wgts <- .Call(C_robustness_weights,
                     x$estimates,
                     length(x$train_ind),
                     rho_cc)

    for (i in seq_along(glbl_wgts)) {
      corrs <- cor(cv_wgts[[i]], glbl_wgts[[i]][x$train_ind, ])
      best_match <- apply(corrs, 2, which.max)

      for (j in seq_along(best_match)) {
        # Preserve the similarity
        matches[[i]][[j]]$similarity[[cv_ind]] <- corrs[best_match[[j]], j]
        # Preserve the prediction error
        matches[[i]][[j]]$predictions[x$test_ind, cvr] <-
          x$estimates[[i]][[best_match[[j]]]]$test_residuals
        # Accumulate the average CV robustness weight
        matches[[i]][[j]]$cv_avg_wgt[x$train_ind, cvr] <-
          matches[[i]][[j]]$cv_avg_wgt[x$train_ind, cvr] +
          cv_wgts[[i]][, best_match[[j]]] / (cv_k - 1L)
      }
    }
  }

  # Summarize prediction errors from matched solutions
  nsol <- sum(lengths(global_ests))
  full_details <- list2DF(list(
    lambda_index   = rep.int(seq_along(global_ests), lengths(global_ests)),
    solution_index = unlist(lapply(lengths(global_ests), seq_len), FALSE, FALSE),
    avg_wrmspe     = numeric(nsol),
    sd_wrmspe      = numeric(nsol),
    avg_wmape      = numeric(nsol),
    sd_wmape       = numeric(nsol),
    avg_wrmspe_cv  = numeric(nsol),
    sd_wrmspe_cv   = numeric(nsol),
    avg_wmape_cv   = numeric(nsol),
    sd_wmape_cv    = numeric(nsol),
    avg_tau_size   = numeric(nsol),
    sd_tau_size    = numeric(nsol),
    avg_similarity = numeric(nsol),
    similarity     = vector('list', nsol)
  ))

  for (rn in seq_len(nsol)) {
    i <- full_details$lambda_index[[rn]]
    j <- full_details$solution_index[[rn]]
    m <- matches[[i]][[j]]

    sum_wgts <- colSums(m$cv_avg_wgt)

    wmape <- colSums(abs(m$predictions) * glbl_wgts[[i]][, j]) / sum(glbl_wgts[[i]][, j])
    wrmspe <- sqrt(colSums(m$predictions^2 * glbl_wgts[[i]][, j]) / sum(glbl_wgts[[i]][, j]))
    wmape_cv <- colSums(abs(m$predictions) * m$cv_avg_wgt) / sum_wgts
    wrmspe_cv <- sqrt(colSums(m$predictions^2 * m$cv_avg_wgt) / sum_wgts)
    tau_size <- apply(m$predictions, 2, tau_size)

    full_details$avg_similarity[[rn]] <- mean(m$similarity)
    full_details$similarity[[rn]]     <- m$similarity
    full_details$avg_tau_size[[rn]]   <- mean(tau_size)
    full_details$sd_tau_size[[rn]]    <- .sd0(tau_size)
    full_details$avg_wrmspe[[rn]]     <- mean(wrmspe)
    full_details$sd_wrmspe[[rn]]      <- .sd0(wrmspe)
    full_details$avg_wrmspe_cv[[rn]]  <- mean(wrmspe_cv)
    full_details$sd_wrmspe_cv[[rn]]   <- .sd0(wrmspe_cv)
    full_details$avg_wmape[[rn]]      <- mean(wmape)
    full_details$sd_wmape[[rn]]       <- .sd0(wmape)
    full_details$avg_wmape_cv[[rn]]   <- mean(wmape_cv)
    full_details$sd_wmape_cv[[rn]]    <- .sd0(wmape_cv)
  }

  # Summarize and select only the "best" solution for each penalization value.
  best_match_global <- list2DF(list(
    lambda_index   = seq_along(global_ests),
    solution_index = integer(length(global_ests)),
    cvavg          = numeric(length(global_ests)),
    cvse           = numeric(length(global_ests))))

  for (i in seq_along(global_ests)) {
    rows <- which(full_details$lambda_index == i)
    best_row <- rows[[which.min(full_details$avg_wrmspe[rows])]]

    best_match_global$solution_index[[i]] <- full_details$solution_index[[best_row]]
    best_match_global$cvavg[[i]]          <- full_details$avg_wrmspe[[best_row]]
    best_match_global$cvse[[i]]           <- full_details$sd_wrmspe[[best_row]]
  }

  list(best = best_match_global, all = full_details)
}

#' @importFrom stats var
.sd0 <- function (x) {
  if (length(x) < 2L) {
    0
  } else {
    sqrt(var(x))
  }
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
    min(candidates[-seq_len(max(which(diff(candidates) > 1)))])
  } else {
    min(candidates)
  }
  type[[best]] <- 'min'
  type[[best_1se]] <- 'se_fact'
  return(type)
}
