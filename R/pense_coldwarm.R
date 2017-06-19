## PENSE with a cold start at the beginning and subsequently uses the previous
## parameter estimate as warm-start
##
## @param lambda_grid a grid of lambda values NOT ADJUSTED for sample size!
## @param start_0 should an initial estimator be computed at the first lambda
##                value or simply initialized at the 0 vector?
#' @importFrom stats mad median
pense_coldwarm <- function(X, y, alpha, lambda_grid, start_0 = FALSE,
                           standardize, pense_options, initest_options,
                           en_options) {
    dX <- dim(X)

    std_data <- standardize_data(X, y, standardize)

    final_estimates <- vector("list", length(lambda_grid))

    ## For the first value of lambda, we will do a cold start
    init.current <- if (isTRUE(start_0)) {
        matrix(c(median(std_data$yc), numeric(dX[2L])), ncol = 1L)
    } else {
        initest_cold(
            std_data$xs,
            std_data$yc,
            alpha,
            lambda_grid[1L],
            pense_options,
            initest_options,
            en_options
        )$initCoef
    }

    for (i in seq_along(lambda_grid)) {
        full_all <- apply(init.current, 2, function(coef, lambda) {
            pen_s_reg(
                std_data$xs,
                std_data$yc,
                alpha,
                lambda,
                init_coef = coef,
                options = pense_options,
                en_options = en_options
            )
        }, lambda = lambda_grid[i])

        best_est <- which.min(sapply(full_all, function(full) {
            full$objF
        }))
        init.current <- matrix(
            c(full_all[[best_est]]$intercept, full_all[[best_est]]$beta),
            ncol = 1L
        )

        final_estimates[[i]] <- full_all[[best_est]]

        if (standardize) {
            final_estimates[[i]]$beta <- final_estimates[[i]]$beta /
                std_data$scale_x
            final_estimates[[i]]$intercept <- final_estimates[[i]]$intercept +
                std_data$muy - drop(final_estimates[[i]]$beta %*% std_data$mux)
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
            init_coef = coef,
            warn = FALSE,
            options = pense_options,
            en_options = en_options
        )

        return(c(
            conc$objF,
            conc$intercept,
            conc$beta
        ))
    })

    objFRanking <- orderOmitTies(initconc[1L, ], tol = initest_options$eps)

    nkeep <- min(length(objFRanking$index), initest_options$keepSolutions)
    indkeep <- objFRanking$index[seq_len(nkeep)]

    return(list(
        initCoef = initconc[-1L, indkeep, drop = FALSE],
        objF = initconc[1L, indkeep]
    ))
}


