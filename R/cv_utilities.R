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

      if (call_with_errors) {
        apply(ordered_predictions - std_data$y, 2, metric)
      } else {
        apply(ordered_predictions, 2, metric, std_data$y)
      }
    })
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
#' @param ncores number of threads to use to match solutions.
#' @param par_cluster parallel cluster to parallelize computations.
#' @param handler_args additional arguments to the handler function.
#' @importFrom Matrix drop
#' @importFrom rlang abort
#' @keywords internal
.run_replicated_cv_ris <- function (std_data, cv_k, cv_repl, cv_est_fun,
                                    global_ests,
                                    min_similarity = 0,
                                    par_cluster = NULL,
                                    rho_cc,
                                    ncores,
                                    handler_args = list()) {
  est_fun <- match.fun(cv_est_fun)
  if (!isTRUE(min_similarity >= 0) || !isTRUE(min_similarity <= 1)) {
    abort("`min_similarity` must be a scalar between 0 and 1.")
  }

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

  dispatcher <- function (test_ind, est_fun, std_data, handler_args) {
    train_x <- std_data$x[-test_ind, , drop = FALSE]
    train_y <- std_data$y[-test_ind]
    test_x <- std_data$x[test_ind, , drop = FALSE]

    train_std <- std_data$cv_standardize(train_x, train_y)

    ## Starting points need to be scaled too
    if (!is.null(handler_args$args$optional_args$individual_starts)) {
      handler_args$args$optional_args$individual_starts <- lapply(
        handler_args$args$optional_args$individual_starts,
        function (starts) { lapply(starts, train_std$standardize_coefs) })
    }
    if (!is.null(handler_args$args$optional_args$shared_starts)) {
      handler_args$args$optional_args$shared_starts <- lapply(
        handler_args$args$optional_args$shared_starts,
        function (starts) { lapply(starts, train_std$standardize_coefs) })
    }

    cv_ests <- est_fun(train_std, test_ind, handler_args)

    train_ind <- seq_along(std_data$y)[-test_ind]

    cv_ests <- lapply(cv_ests, function (ests_lambda) {
      lapply(ests_lambda, function (est) {
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
  # cv_ests <- cl_handler(test_segments,
  #                       dispatcher,
  #                       est_fun = est_fun,
  #                       std_data = std_data,
  #                       handler_args = handler_args)
  .RISCV_DEBUG_CACHE_PATH <- Sys.getenv("RISCV_DEBUG_CACHE", "")

  if (nzchar(.RISCV_DEBUG_CACHE_PATH) && file.exists(.RISCV_DEBUG_CACHE_PATH)) {
    cv_ests <- readRDS(.RISCV_DEBUG_CACHE_PATH)
  } else {
    cv_ests <- cl_handler(test_segments,
                          dispatcher,
                          est_fun = est_fun,
                          std_data = std_data,
                          handler_args = handler_args)

    if (nzchar(.RISCV_DEBUG_CACHE_PATH)) {
      saveRDS(cv_ests, file = .RISCV_DEBUG_CACHE_PATH)
    }
  }

  matches <- .Call(C_match_solutions_by_weight,
                   cv_ests,
                   global_ests,
                   as.integer(cv_k),
                   length(std_data$y),
                   rho_cc,
                   ncores)

  best_match_global <- as.data.frame(t(
    vapply(matches, FUN.VALUE = numeric(4), FUN = function (lambda_match) {
      avgs <- vapply(lambda_match, FUN.VALUE = numeric(3), FUN = function (sol_match) {
        c(cvavg = mean(sol_match$wmspe),
          cvse = if (length(sol_match$wmspe) > 1L) { sd(sol_match$wmspe) } else { 0 },
          avg_similarity = median(sol_match$rankcorr))
      })
      avgs <- rbind(avgs, solution_index = seq_len(ncol(avgs)))

      candidates <- which(avgs['avg_similarity', ] >= min_similarity)
      if (length(candidates) > 0L) {
        avgs <- avgs[, candidates, drop = FALSE]
      }

      avgs[, which.min(avgs['cvavg', ])]
    })))
  best_match_global$lambda_index <- seq_len(nrow(best_match_global))

  full_details <- do.call(rbind, lapply(seq_along(matches), function (lambda_ind) {
    lambda_match <- matches[[lambda_ind]]

    list2DF(list(lambda_index = rep.int(lambda_ind, length(lambda_match)),
                 solution_index = seq_along(lambda_match),
                 avg_wmspe = vapply(lambda_match, FUN.VALUE = numeric(1),
                                    FUN = function (sol_match) {
                                      mean(sol_match$wmspe)
                                    }),
                 sd_wmspe = vapply(lambda_match, FUN.VALUE = numeric(1),
                                   FUN = function (sol_match) {
                                     if (length(sol_match$wmspe) > 1L) {
                                       sd(sol_match$wmspe)
                                     } else {
                                       0
                                     }
                                   }),
                 avg_similarity = vapply(lambda_match, FUN.VALUE = numeric(1),
                                         FUN = function (sol_match) {
                                           median(sol_match$rankcorr)
                                         }),
                 rankcorr = lapply(lambda_match, `[[`, 'rankcorr')))
  }))

  list(best = best_match_global, all = full_details)
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
