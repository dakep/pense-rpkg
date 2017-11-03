## Get (one or more) cold initial estimates at the given lambda.
## @param lambda penalty parameter
#' @importFrom stats mad median
#' @importFrom Matrix Matrix sparseMatrix drop
#' @importClassesFrom Matrix dgCMatrix
pense_init_cold <- function(x, y, alpha, lambda, standardize,
                            pense_options, initest_options,
                            en_options) {
    dx <- dim(x)
    std_data <- standardize_data(x, y, standardize)
    lambda <- lambda / max(std_data$scale_x)

    # Compute cold initial estimate
    initraw <- enpy(
        std_data$xs,
        std_data$yc,
        alpha,
        lambda,
        delta = pense_options$bdp,
        cc = pense_options$cc,
        options = initest_options,
        en_options = en_options
    )

    # Only do a few refinement steps
    pense_options$maxit <- initest_options$maxitPenseRefinement

    ## Compute a "short" PENSE for each candidate solution
    initconc <- apply(initraw$coeff, 2, function(coef) {
        conc <- pen_s_reg(
            std_data$xs,
            std_data$yc,
            alpha,
            lambda,
            init_int = coef[1L],
            init_coef = coef[-1L],
            warn = FALSE,
            options = pense_options,
            en_options = en_options
        )

        return(conc[c("intercept", "beta", "objF")])
    })

    # Select only the best solutions
    obj_vals <- vapply(
        initconc,
        "[[",
        "objF",
        FUN.VALUE = numeric(1L),
        USE.NAMES = FALSE
    )
    objFRanking <- orderOmitTies(obj_vals, tol = initest_options$eps)

    nkeep <- min(length(objFRanking$index), initest_options$keepSolutions)
    indkeep <- objFRanking$index[seq_len(nkeep)]

    return(lapply(initconc[indkeep], std_data$unstandardize_coefs))
}


## Compute the fully iterated PENSE solution using the solution at the previous
## lambda value as well as any given initial estimates in `initial_ests`.
##
## @param lambda_grid a grid of lambda values NOT ADJUSTED for sample size!
## @param initial_ests a list the same length as lambda_grid. If the first
##      item is NULL, it will be replaced by a 0-start (the behavior is the same
##      when the argument is missing)
## @param refine_it number of iterations to refine the initial estimates
## @param nkeep number of initial estimates kept after refining
#' @importFrom stats mad median
#' @importFrom Matrix Matrix sparseMatrix drop
#' @importClassesFrom Matrix dgCMatrix
pense_full <- function(x, y, alpha, lambda_grid, standardize,
                       initial_ests, refine_it, nkeep,
                       pense_options, en_options, warn = TRUE) {
    dx <- dim(x)
    std_data <- standardize_data(x, y, standardize)
    lambda_grid <- lambda_grid / max(std_data$scale_x)
    final_estimates <- vector("list", length(lambda_grid))

    if (missing(initial_ests)) {
        initial_ests <- vector("list", length(lambda_grid))
    }

    init_current <- vector("list", 0L)

    ## If the first initial_ests is NULL or empty, use a 0-start
    if (is.null(initial_ests[[1L]]) || length(initial_ests[[1L]]) == 0) {
        init_current <- list(
            list(
                intercept = median(std_data$yc),
                beta = sparseMatrix(
                    i = integer(0L),
                    j = integer(0L),
                    x = numeric(0L),
                    dims = c(dx[2L], 1L)
                )
            )
        )
    }

    # After 5 iterations, we only keep the top 5 estimators
    original_maxit <- pense_options$maxit
    refine_it <- as.integer(refine_it)
    refine_nkeep <- as.integer(nkeep)

    for (i in seq_along(lambda_grid)) {
        init_current <- c(init_current, initial_ests[[i]])
        pense_options$maxit <- refine_it

        half_all <- lapply(init_current, function (ic) {
            pen_s_reg(
                std_data$xs,
                std_data$yc,
                alpha = alpha,
                lambda = lambda_grid[[i]],
                init_int = ic$intercept,
                init_coef = ic$beta,
                options = pense_options,
                en_options = en_options,
                warn = FALSE
            )
        })

        half_objf <- sapply(half_all, '[[', 'objF')
        objf_ranking <- orderOmitTies(half_objf, pense_options$eps)
        nkeep <- min(length(objf_ranking$index), refine_nkeep)
        indkeep <- objf_ranking$index[seq_len(nkeep)]

        pense_options$maxit <- original_maxit

        full_kept <- lapply(half_all[indkeep], function (ic) {
            pen_s_reg(
                std_data$xs,
                std_data$yc,
                alpha = alpha,
                lambda = lambda_grid[[i]],
                init_int = ic$intercept,
                init_coef = ic$beta,
                options = pense_options,
                en_options = en_options,
                warn = FALSE
            )
        })

        best_est <- which.min(sapply(full_kept, '[[', 'objF'))
        init_current <- list(full_kept[[best_est]][c("intercept", "beta")])

        if (!isTRUE(full_kept[[best_est]]$rel_change < pense_options$eps) && isTRUE(warn)) {
            warning(sprintf("PENSE did not converge for lambda = %.6f", lambda_grid[[i]]))
        }

        final_estimates[[i]] <- std_data$unstandardize_coefs(full_kept[[best_est]])
    }

    return(final_estimates)
}


initest_cold <- function(x, y, alpha, lambda, pense_options,
                         initest_options, en_options) {
    initraw <- enpy(
        x,
        y,
        alpha,
        lambda,
        delta = pense_options$bdp,
        cc = pense_options$cc,
        options = initest_options,
        en_options = en_options
    )

    short_pense_options <- pense_options
    short_pense_options$maxit <- 10L

    ## Compute a "short" PENSE for each candidate solution
    initconc <- apply(initraw$coeff, 2, function(coef) {
        conc <- pen_s_reg(
            x,
            y,
            alpha,
            lambda,
            init_int = coef[1L],
            init_coef = coef[-1L],
            warn = FALSE,
            options = short_pense_options,
            en_options = en_options
        )

        return(conc[c("intercept", "beta", "objF")])
    })

    obj_vals <- vapply(initconc, "[[", "objF", FUN.VALUE = numeric(1L),
                       USE.NAMES = FALSE)
    objFRanking <- orderOmitTies(obj_vals, tol = initest_options$eps)

    nkeep <- min(length(objFRanking$index), initest_options$keepSolutions)
    indkeep <- objFRanking$index[seq_len(nkeep)]

    return(initconc[indkeep])
}


