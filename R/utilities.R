nameCoefVec <- function(coef, X) {
    dn <- dimnames(X)
    xnames <- paste("X", seq_len(ncol(X)), sep = "")

    if (!is.null(dn) && !is.null(dn[[2L]])) {
        xnames <- dn[[2L]]
    }

    if (is.matrix(coef)) {
        rownames(coef) <- c("(Intercept)", xnames)
    } else {
        names(coef) <- c("(Intercept)", xnames)
    }
    return(coef)
}

orderOmitTies <- function(x, tol) {
    ord.x <- sort.list(x, na.last = NA, method = "quick")
    sorted <- x[ord.x]
    diffs <- diff(sorted)

    filtered.x <- c(sorted[1L], sorted[-1L][diffs > tol])
    filtered.ind <- c(ord.x[1L], ord.x[-1L][diffs > tol])

    return(list(
        clean = filtered.x,
        index = filtered.ind
    ))
}

## Get the constant needed for consistency for the given delta
## and the given rho function
#' @importFrom robustbase .Mchi
consistency.rho <- function(delta, int.rho.fun, interval = c(0.3, 10)) {
    if (is.character(int.rho.fun)) {
        int.rho.fun <- switch (int.rho.fun,
                               huber = 0L,
                               bisquare = 1L,
                               1L
        )
    }

    integrand <- function(x, cc) {
        dnorm(x) * .Mchi(x, cc, int.rho.fun)
    }

    expectation <- if (int.rho.fun == 1L) {
        # For bisquare we have the closed form solution to the expectation as
        function(cc, delta) {
            pnorm.mcc <- 2 * pnorm(-cc)
            1/cc^6 * exp(-(cc^2/2)) * (
                -cc * (15 - 4 * cc^2 + cc^4) * sqrt(2 / pi) +
                    3 * (5 - 3 * cc^2 + cc^4) * exp(cc^2/2) * (1 - pnorm.mcc) +
                    cc^6 * exp(cc^2/2) * pnorm.mcc
            ) - delta
        }
    } else {
        function(cc, delta) {
            integrate(integrand, lower = -Inf, upper = Inf, cc)$value - delta
        }
    }

    uniroot(expectation, interval = interval, delta)$root
}


standardize_data <- function (x, y, standardize) {
    ret_list <- list(
        scale_x = 1,
        mux = 0,
        muy = 0,
        xs = x,
        yc = y
    )

    ## standardize data
    if (isTRUE(standardize)) {
        ret_list$scale_x <- apply(x, 2, mad)
        ret_list$mux <- apply(x, 2, median)
        ret_list$muy <- median(y)

        ret_list$xs <- scale(x, center = ret_list$mux, scale = ret_list$scale_x)
        ret_list$yc <- y - ret_list$muy
    }

    return(ret_list)
}

