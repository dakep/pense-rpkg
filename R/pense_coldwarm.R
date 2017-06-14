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

    scale.x <- 1
    mux <- 0
    muy <- 0

    if (isTRUE(standardize)) {
        scale.x <- apply(X, 2, mad)
        mux <- apply(X, 2, median)
        muy <- median(y)

        X <- scale(X, center = mux, scale = scale.x)
        y <- y - muy
    }

    final.estimates <- vector("list", length(lambda_grid))

    ## For the first value of lambda, we will do a cold start
    init.current <- if (isTRUE(start_0)) {
        matrix(c(median(y), numeric(dX[2L])), ncol = 1L)
    } else {
        initest_cold(X, y, alpha, lambda_grid[1L], pense_options,
                     initest_options, en_options)$initCoef
    }

    for (i in seq_along(lambda_grid)) {
        full.all <- apply(init.current, 2, function(coef, lambda) {
            pen_s_reg(
                X,
                y,
                alpha,
                lambda,
                init_coef = coef,
                options = pense_options,
                en_options = en_options
            )
        }, lambda = lambda_grid[i])

        best.est <- which.min(sapply(full.all, function(full) {
            full$objF
        }))
        init.current <- matrix(
            c(full.all[[best.est]]$intercept, full.all[[best.est]]$beta),
            ncol = 1L
        )

        final.estimates[[i]] <- full.all[[best.est]]

        if (standardize == TRUE) {
            final.estimates[[i]]$beta <- final.estimates[[i]]$beta / scale.x
            final.estimates[[i]]$intercept <- final.estimates[[i]]$intercept +
                muy - drop(final.estimates[[i]]$beta %*% mux)
        }
    }

    return(final.estimates)
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


