#' @importFrom lifecycle deprecated is_present deprecate_stop deprecate_stop
#' @importFrom rlang warn abort is_missing
#' @importFrom methods is
#' @importFrom stats runif
.pense_args <- function (...,
                         options = deprecated(),
                         init_options = deprecated(),
                         en_options = deprecated(),
                         initial = deprecated(),
                         warm_reset = deprecated()) {
  ## Check presence of deprecated arguments
  # Translate options for initial estimates
  if (is_present(warm_reset)) {
    deprecate_stop('2.0.0', 'pense(warm_reset=)', 'pense(nlambda_enpy=)')
  }
  if (is_present(initial)) {
    deprecate_stop('2.0.0', 'pense(initial=)',
                   details = paste("Use arguments `nlambda_enpy`,",
                                   "`enpy_specific`, `add_zero_based`",
                                   "and `other_starts` instead."))
  }

  # Translate algorithm options
  if (is_present(en_options)) {
    deprecate_stop('2.0.0', 'pense(en_options=)',
                   details = paste("The LS-EN algorithm is specified",
                                   "with arguments `algorithm_opts`",
                                   "and `enpy_opts`."))
  }
  if (is_present(init_options)) {
    deprecate_stop('2.0.0', 'pense(init_options=)', 'pense(enpy_opts=)')
  }
  if (is_present(options)) {
    deprecate_stop('2.0.0', 'pense(options=)',
                   details = paste("Use arguments `bdp`, `eps`,",
                                   "`mscale_opts` and `algorithm_opts`",
                                   "to `pense()` instead."))
  }

  ## Pull in default values for the arguments from the pense() function
  ## definition.
  defaults <- formals(pense)
  args <- list2env(list(...))
  for (argname in names(defaults)) {
    if (!(exists(argname, args, inherits = FALSE))) {
      # evaluate the default values in this environment
      args[[argname]] <- if (!is_missing(defaults[[argname]])) {
        eval(defaults[[argname]], envir = args)
      } else {
        defaults[[argname]]
      }
    }
  }

  args <- as.list(args)
  args$optional_args <- list()

  ## Process input arguments
  response <- .validate_response(args$y)
  args$y <- response$values
  args$binary_response <- response$binary
  x_dim <- dim(args$x)

  if (length(args$y) != x_dim[[1L]]) {
    abort("Number of observations in `x` and `y` does not match.")
  } else if (x_dim[[2L]] <= 1L) {
    abort("`x` must be a matrix with at least 2 columns.")
  }

  args$alpha <- .as(args$alpha, 'numeric')
  if (any(args$alpha < 0 | args$alpha > 1)) {
    abort("`alpha` is outside 0 and 1.")
  } else if (any(args$alpha < sqrt(.Machine$double.eps))) {
    args$alpha[which(args$alpha < sqrt(.Machine$double.eps))] <- 0
    if (any(args$alpha > 0)) {
      abort("`alpha=0` cannot be mixed with other `alpha` values.")
    }
  }

  if (!is_missing(args$enpy_lambda)) {
    args$nlambda_enpy <- length(args$enpy_lambda)
  }

  args$pense_opts <- with(
    args,
    list(algo_opts = algorithm_opts,
         strategy_0 = isTRUE(add_zero_based),
         strategy_enpy_individual = isTRUE(enpy_specific) && (nlambda_enpy > 0L),
         strategy_enpy_shared = !isTRUE(enpy_specific) && (nlambda_enpy > 0L),
         strategy_other_individual = FALSE,
         strategy_other_shared = FALSE,
         algorithm = .pense_algorithm_id(algorithm_opts),
         intercept = !isFALSE(intercept),
         warm_starts = !isFALSE(carry_forward[[1L]]),
         eps = .as(eps[[1L]], 'numeric'),
         comparison_tol = .as(comparison_tol[[1L]], 'numeric'),
         explore_tol = .as(explore_tol[[1L]], 'numeric'),
         explore_it = .as(explore_it[[1L]], 'integer'),
         nr_tracks = .as(explore_solutions[[1L]], 'integer'),
         max_optima = .as(max_solutions[[1L]], 'integer'),
         num_threads = max(1L, .as(ncores[[1L]], 'integer')),
         sparse = isTRUE(sparse),
         mscale = .full_mscale_algo_options(bdp = bdp, cc = cc,
                                            mscale_opts = mscale_opts)))

  if (args$pense_opts$explore_tol < args$pense_opts$eps) {
    abort("`explore_tol` must not be less than `eps`")
  }
  if (args$pense_opts$comparison_tol < args$pense_opts$eps) {
    abort("`comparison_tol` must not be less than `eps`")
  }
  if (args$pense_opts$explore_it < 0L) {
    abort("`explore_it` must not be less than 0")
  }

  # Check EN algorithm for ENPY
  args$enpy_opts$en_options <- .select_en_algorithm(args$enpy_opts$en_options,
                                                    args$alpha,
                                                    args$pense_opts$sparse,
                                                    args$eps)
  args$pense_opts$sparse <- args$enpy_opts$en_options$sparse

  # If using the MM algorithm, ensure that the EN options are set.
  if (identical(args$pense_opts$algorithm, .k_pense_algo_mm)) {
    args$pense_opts$algo_opts$en_options <-
      .select_en_algorithm(args$pense_opts$algo_opts$en_options,
                           args$alpha,
                           args$pense_opts$sparse,
                           args$eps)
    if (!isTRUE(args$pense_opts$sparse ==
                args$pense_opts$algo_opts$en_options$sparse)) {
      abort(paste("The `sparse` option for the EN-PY algorithm and the",
                  "MM algorithm for PENSE disagree."))
    }
  }

  # Set the number of cores for the ENPY options
  if (args$pense_opts$num_threads > 1L && !isTRUE(.k_multithreading_support)) {
    warn("Multithreading not supported. Using only 1 core.")
    args$pense_opts$num_threads <- 1L
  }
  args$enpy_opts$num_threads <- args$pense_opts$num_threads

  # Standardizing the data
  standardize <- if (is.character(args$standardize)) {
    if (pmatch(args$standardize[[1L]], 'cv_only', nomatch = 0L) == 1L) {
      args$standardize <- 'cv_only'
    } else {
      abort("`standardize` must be either TRUE/FALSE or \"cv_only\".")
    }
  } else {
    isTRUE(args$standardize)
  }

  # Check penalty loadings
  if (!is_missing(args$penalty_loadings) && !is.null(args$penalty_loadings)) {
    checked_pls <- .prepare_penalty_loadings(args$penalty_loadings,
                                             x = args$x,
                                             alpha = args$alpha,
                                             sparse = args$pense_opts$sparse)
    args$penalty_loadings <- checked_pls$loadings
    args$restore_coef_length <- checked_pls$restore_fun
    args$x <- checked_pls$trimmed_x
  } else {
    args$restore_coef_length <- function (coef) coef
    args$penalty_loadings <- NULL
  }

  if (ncol(args$x) == 0L) {
    args$pense_opts$intercept <- TRUE
    warn(paste("All values in `penalty_loadings` are infinite.",
               "Only computing the intercept."))

    args$std_data <- .standardize_data(
      matrix(runif(x_dim[[1L]]), ncol = 1L),
      args$y,
      intercept = TRUE,
      sparse = args$pense_opts$sparse,
      standardize = standardize,
      robust = TRUE,
      mscale_opts = args$mscale_opts,
      bdp = args$pense_opts$mscale$delta,
      cc = args$pense_opts$mscale$cc)

    # Compute only the 0-based solution.
    args$pense_opts$strategy_enpy_individual <- FALSE
    args$pense_opts$strategy_enpy_shared <- FALSE
    args$pense_opts$strategy_0 <- TRUE
    args$lambda <- lapply(args$alpha,
                          FUN = .pense_lambda_grid,
                          x = args$std_data$x,
                          y = args$std_data$y,
                          nlambda = 1,
                          lambda_min_ratio = 1,
                          pense_options = args$pense_opts,
                          penalty_loadings = NULL)
    args$enpy_lambda_inds <- rep(list(integer(0L)), length(args$alpha))

    return(args)
  }

  args$std_data <- .standardize_data(
    args$x,
    args$y,
    intercept = args$pense_opts$intercept,
    standardize = standardize,
    robust = TRUE,
    sparse = args$pense_opts$sparse,
    mscale_opts = args$mscale_opts,
    bdp = args$pense_opts$mscale$delta,
    cc = args$pense_opts$mscale$cc)

  # Scale penalty loadings appropriately
  args$penalty_loadings <- args$penalty_loadings / args$std_data$scale_x
  if (length(args$penalty_loadings) == 0L) {
    args$penalty_loadings <- NULL
  }

  # Determine lambda grid
  args$lambda <- if (is_missing(args$lambda) || is.null(args$lambda)) {
    if (is_missing(args$lambda_min_ratio)) {
      args$lambda_min_ratio <- NULL
    }

    lapply(args$alpha,
           FUN = .pense_lambda_grid,
           x = args$std_data$x,
           y = args$std_data$y,
           nlambda = args$nlambda,
           lambda_min_ratio = args$lambda_min_ratio,
           pense_options = args$pense_opts,
           penalty_loadings = args$penalty_loadings)
  } else if (!is.list(args$lambda)) {
    rep.int(list(sort(.as(args$lambda, 'numeric'), decreasing = TRUE)),
            length(args$alpha))
  } else if (identical(length(args$lambda), length(args$alpha))) {
    lapply(args$lambda, function (l) {
      sort(.as(l, 'numeric'), decreasing = TRUE)
    })
  } else {
    abort(paste("`lambda` must either be a numeric vector or a list the",
                "same length as `alpha`."))
  }

  # Split the `other_starts` into individual and shared starts.
  if (!is_missing(args$other_starts)) {
    if (is(args$other_starts, 'starting_point')) {
      args$other_starts <- structure(list(args$other_starts),
                                     class = 'starting_points')
    } else if (!is(args$other_starts, 'starting_points')) {
      abort(paste("`other_starts` must be a list of starting points created by",
                  "`starting_point()`, `enpy_initial_estimates()`,",
                  "or a combination thereof."))
    }

    # Identify which other starts are shared and which are specific.
    other_starts_shared <- vapply(args$other_starts,
                                  FUN.VALUE = logical(1L),
                                  FUN = is,
                                  'shared_starting_point')
    other_starts_specific <- vapply(args$other_starts,
                                    FUN.VALUE = logical(1L),
                                    FUN = is,
                                    'specific_starting_point')

    # Ensure the `beta` coefficients in `other_starts` agree with the
    # desired vector class (sparse vs. dense) and standardize them.
    other_starts <- lapply(
      .sparsify_other_starts(args$other_starts, args$pense_opts$sparse),
      args$std_data$standardize_coefs)

    if (any(other_starts_shared)) {
      args$pense_opts$strategy_other_shared <- TRUE
      args$optional_args$shared_starts <- other_starts[other_starts_shared]
    }
    if (any(other_starts_specific)) {
      args$pense_opts$strategy_other_individual <- TRUE
      ind_starts <- list()
      for (ai in seq_along(args$alpha)) {
        new_ind_starts <- .make_initest_list(
          other_starts[other_starts_specific],
          lambda = args$lambda[[ai]],
          alpha = args$alpha[[ai]],
          sparse = args$pense_opts$sparse)

        ind_starts <- c(ind_starts, new_ind_starts$starting_points)
        args$lambda[[ai]] <- new_ind_starts$extended_lambda
      }
      args$optional_args$individual_starts <- ind_starts
    }
  }

  # Determine ENPY lambda grid
  args$enpy_lambda_inds <- if (args$pense_opts$strategy_enpy_individual ||
                               args$pense_opts$strategy_enpy_shared) {
    if (is_missing(args$enpy_lambda) || is.null(args$enpy_lambda)) {
      lapply(args$lambda, function (l) {
        nlambda_enpy <- min(length(l), args$nlambda_enpy)
        eq_spaced_seq <- seq(1, length(l), length.out = nlambda_enpy + 1L)
        as.integer(ceiling(eq_spaced_seq)[-(nlambda_enpy + 1L)])
      })
    } else if (is.list(args$enpy_lambda)) {
      mapply(args$enpy_lambda, args$lambda,
             FUN = .approx_match,
             SIMPLIFY = FALSE, USE.NAMES = FALSE)
    } else if (is.numeric(args$enpy_lambda)) {
      vapply(args$lambda,
             FUN.VALUE = integer(length(args$enpy_lambda)),
             FUN = .approx_match,
             x = args$enpy_lambda)
    } else {
      abort(paste("`enpy_lambda` must either be a numeric vector or",
                  "a list the same length as `alpha`."))
    }
  } else {
    rep(list(integer(0L)), length(args$alpha))
  }

  # Extend lambda grids if necessary
  for (ai in seq_along(args$alpha)) {
    if (anyNA(args$enpy_lambda_inds[[ai]])) {
      if (is.list(args$enpy_lambda)) {
        args$lambda[[ai]] <- sort(c(args$lambda[[ai]], args$enpy_lambda[[ai]]),
                                  decreasing = TRUE)
        args$enpy_lambda_inds[[ai]] <- .approx_match(args$enpy_lambda[[ai]],
                                                     args$lambda[[ai]])
      } else if (is.numeric(args$enpy_lambda)) {
        args$lambda[[ai]] <- sort(c(args$lambda[[ai]], args$enpy_lambda),
                                  decreasing = TRUE)
        args$enpy_lambda_inds[[ai]] <- .approx_match(args$enpy_lambda,
                                                     args$lambda[[ai]])
      }
    }
  }

  for (lambda_grid in args$lambda) {
    if (any(lambda_grid < .Machine$double.eps)) {
      abort("All values in `lambda` must be positive.")
    }
  }

  args
}

## Make a list of initial estimates
##
## @return a list with 2 components:
##   `extended_lambda` an extended grid of penalization levels to contain
##                     both the given lambda values plus the lambda values
##                     in `other_starts`.
##   `starting_points` a list the same length as `extended_lambda` with a list
##                     of initial estimates for each value in `extended_lambda`.
#' @importFrom rlang warn
.make_initest_list <- function (other_starts, lambda, alpha, sparse) {
  if (length(other_starts) == 0L) {
    return(list(extended_lambda = lambda,
                starting_points = rep.int(list(list()), length(lambda))))
  }

  # Check for wrong starting points without lambda
  init_est_lambda <- unlist(lapply(other_starts, `[[`, 'lambda'),
                            use.names = FALSE, recursive = FALSE)

  if (length(init_est_lambda) != length(other_starts)) {
    abort(paste("Some starting points in `other_starts` are marked as",
                "\"specific\" but do not have a `lambda` component."))
  }

  init_est_lambda <- .as(init_est_lambda, 'numeric')

  # Check for wrong starting points without alpha
  init_est_alpha <- unlist(lapply(other_starts, `[[`, 'alpha'),
                           use.names = FALSE, recursive = FALSE)

  if (length(init_est_alpha) != length(other_starts)) {
    abort(paste("Some starting points in `other_starts` are marked as",
                "\"specific\" but do not have an `alpha` component."))
  }

  init_est_alpha <- .as(init_est_alpha, 'numeric')
  correct_alpha <- which(abs(init_est_alpha - alpha) <
                           sqrt(.Machine$double.eps))

  if (length(correct_alpha) == 0L) {
    return(list(extended_lambda = lambda,
                starting_points = rep.int(list(list()), length(lambda))))
  }
  init_est_lambda <- init_est_lambda[correct_alpha]
  other_starts <- other_starts[correct_alpha]

  init_est_inds <- .approx_match(init_est_lambda, lambda)
  new_initest_lambda <- which(is.na(init_est_inds))

  if (length(new_initest_lambda) > 0L) {
    # Some starts are for unknown lambda. Add lambdas to the grid!
    lambda <- sort(c(lambda, unique(init_est_lambda[new_initest_lambda])),
                   decreasing = TRUE)
    init_est_inds <- .approx_match(init_est_lambda, lambda)
  }

  starting_points <- lapply(seq_along(lambda), function (i) {
    matches <- which(i == init_est_inds)
    if (length(matches) > 0L) {
      return(other_starts[matches])
    } else {
      return(list())
    }
  })

  return(list(extended_lambda = lambda, starting_points = starting_points))
}

## Get the smallest lambda such that the PENSE estimate gives the empty model.
.pense_max_lambda <- function (x, y, alpha, pense_options,
                               penalty_loadings = NULL) {
  optional_args <- list()
  if (!is.null(penalty_loadings)) {
    optional_args$pen_loadings <- penalty_loadings
  }
  .Call(C_pense_max_lambda, x, y, pense_options, optional_args) /
    max(0.01, alpha)
}

## Generate a log-spaced grid of decreasing lambda values
#' @importFrom rlang abort
.pense_lambda_grid <- function (x, y, alpha, nlambda, lambda_min_ratio,
                                pense_options, penalty_loadings) {
  alpha <- max(0.01, alpha)
  x_dim <- dim(x)
  if (is.null(lambda_min_ratio)) {
    lambda_min_ratio <- alpha * if (x_dim[[1L]] > x_dim[[2L]]) {
      1e-3
    } else {
      1e-2
    }
  }
  max_lambda <- .pense_max_lambda(x, y, alpha, pense_options, penalty_loadings)

  if (!isTRUE(max_lambda > .Machine$double.eps)) {
    abort("Cannot determine maximum lambda. Scale of response is likely 0.")
  }

  rev(exp(seq(log(lambda_min_ratio * max_lambda), log(max_lambda),
              length.out = nlambda)))
}
