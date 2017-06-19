## PENSE with a cold start at the beginning and subsequently uses the previous
## parameter estimate as warm-start
##
## @param lambda_grid a grid of lambda values NOT ADJUSTED for sample size!
## @param start_0 should an initial estimator be computed at the first lambda
##                value or simply initialized at the 0 vector?
#' @importFrom stats mad median
#' @importFrom Matrix Matrix sparseMatrix drop
#' @importClassesFrom Matrix dgCMatrix
pense_coldwarm <- function(X, y, alpha, lambda_grid, start_0 = FALSE,
                           standardize, pense_options, initest_options,
                           en_options) {
    dX <- dim(X)

    std_data <- standardize_data(X, y, standardize)

    final_estimates <- vector("list", length(lambda_grid))

    ## For the first value of lambda, we will do a cold start
    if (isTRUE(start_0)) {
        init_current <- list(
            list(
                intercept = median(std_data$yc),
                beta = sparseMatrix(
                    i = integer(0L),
                    j = integer(0L),
                    x = numeric(0L),
                    dims = c(dX[2L], 1L)
                )
            )
        )
    } else {
        init_current <- initest_cold(
            std_data$xs,
            std_data$yc,
            alpha,
            lambda_grid[1L],
            pense_options,
            initest_options,
            en_options
        )
    }

    for (i in seq_along(lambda_grid)) {
        full_all <- lapply(init_current, function (ic) {
            pen_s_reg(
                std_data$xs,
                std_data$yc,
                alpha = alpha,
                lambda = lambda_grid[i],
                init_int = ic$intercept,
                init_coef = ic$beta,
                options = pense_options,
                en_options = en_options
            )
        })

        best_est <- which.min(sapply(full_all, function(full) {
            full$objF
        }))

        init_current <- list(full_all[[best_est]][c("intercept", "beta")])

        final_estimates[[i]] <- full_all[[best_est]]

        if (standardize) {
            final_estimates[[i]]$beta <- final_estimates[[i]]$beta /
                std_data$scale_x
            final_estimates[[i]]$intercept <- final_estimates[[i]]$intercept +
                std_data$muy -
                drop(std_data$mux %*% final_estimates[[i]]$beta)
        }
    }

    return(final_estimates)
}

initest_cold <- function(X, y, alpha, lambda, pense_options,
                         initest_options, en_options) {
    initraw <- enpy(
        X,
        y,
        alpha,
        lambda,
        options = initest_options,
        en_options = en_options
    )

    ## Compute a "short" PENSE for each candidate solution
    initconc <- apply(initraw$coeff, 2, function(coef) {
        conc <- pen_s_reg(
            X,
            y,
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

    obj_vals <- vapply(initconc, "[[", "objF", FUN.VALUE = numeric(1L),
                       USE.NAMES = FALSE)
    objFRanking <- orderOmitTies(obj_vals, tol = initest_options$eps)

    nkeep <- min(length(objFRanking$index), initest_options$keepSolutions)
    indkeep <- objFRanking$index[seq_len(nkeep)]

    return(initconc[indkeep])
}


