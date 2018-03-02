#' @importFrom methods is
nameCoefVec <- function(coef, x) {
    dn <- dimnames(x)
    xnames <- paste("X", seq_len(ncol(x)), sep = "")

    if (!is.null(dn) && !is.null(dn[[2L]])) {
        xnames <- dn[[2L]]
    }

    if (is.matrix(coef) || is(coef, "Matrix")) {
        rownames(coef) <- c("(Intercept)", xnames)
    } else {
        names(coef) <- c("(Intercept)", xnames)
    }
    return(coef)
}

## Get the order of a list of values removing the duplicate
## entries
##
orderOmitTies <- function(x, tol = 1e-6) {
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
#' @importFrom stats dnorm pnorm integrate uniroot
consistency.rho <- function(delta, int.rho.fun) {
    if (is.character(int.rho.fun)) {
        int.rho.fun <- .rho2IntRho(int.rho.fun)
    }

    ##
    ## Pre-computed values for some delta values
    ##
    if (abs(delta - 0.5) < sqrt(.Machine$double.eps)) {
        return(switch(
            as.character(int.rho.fun),
            "0" = 1.3684820, # huber
            "1" = 1.5476450, # bisquare
            "5" = 0.5773503 # gauss
        ))
    } else if (abs(delta - 0.25) < sqrt(.Machine$double.eps)) {
        return(switch(
            as.character(int.rho.fun),
            "0" = 1.988013, # huber
            "1" = 2.937015, # bisquare
            "5" = 1.133893  # gauss
        ))
    } else if (abs(delta - 0.1) < sqrt(.Machine$double.eps)) {
        return(switch(
            as.character(int.rho.fun),
            "0" = 3.161931, # huber
            "1" = 5.182361, # bisquare
            "5" = 2.064742  # gauss
        ))
    } else if (delta < 0.005) {
        return(50) # ~.1% bdp for bisquare, 9.6e-5% for huber, 0.02% for gauss
    }

    integrand_huber <- function(x, cc) {
        dnorm(x) * .Mchi(x, cc, 0L) / (0.5 * cc * cc)
    }
    integrand_gauss <- function(x, cc) {
        dnorm(x) * -expm1(-((x * x) / (cc * cc)) * 0.5)
    }

    if (int.rho.fun == 1L) {
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
    } else if (int.rho.fun == 0L) {
        integral_interval <- if (delta > 0.1) {
            c(.1, 7)
        } else {
            c(3, 30)
        }
        expectation <- function(cc, delta) {
            integrate(integrand_huber, lower = -Inf, upper = Inf, cc)$value - delta
        }
    } else if (int.rho.fun == 5L) {
        integral_interval <- if (delta > 0.1) {
            c(.5, 2.5)
        } else {
            c(2, 10)
        }

        expectation <- function(cc, delta) {
            integrate(integrand_gauss, lower = -Inf, upper = Inf, cc)$value - delta
        }
    }

    uniroot(expectation, interval = integral_interval, delta)$root
}


## Standardize data (depending on the `standardize` parameter)
##
## The returned list has the information on the standardization
## as well as functions to (un)standardize regression coefficients
##
#' @importFrom Matrix drop
standardize_data <- function (x, y, standardize, robust = TRUE) {
    ret_list <- list(
        scale_x = 1,
        mux = 0,
        muy = 0,
        xs = x,
        yc = y
    )

    ## standardize data
    if (isTRUE(standardize)) {
        if (!isTRUE(robust)) {
            ret_list$scale_x <- apply(x, 2, sd)
            ret_list$mux <- colMeans(x)
            ret_list$muy <- mean(y)
        } else {
            ret_list$scale_x <- apply(x, 2, mad)
            ret_list$mux <- apply(x, 2, median)
            ret_list$muy <- median(y)
        }

        if (!isTRUE(all(ret_list$scale_x > 0))) {
            stop("One or more variables in x have a MAD of 0. Can not use ",
                 "`standardize = TRUE`!")
        }

        ret_list$xs <- scale(x, center = ret_list$mux, scale = ret_list$scale_x)
        ret_list$yc <- y - ret_list$muy
    }

    ret_list$standardize_coefs <- function(coef_obj) {
        if (!isTRUE(standardize)) {
            return(coef_obj)
        }

        coef_obj$intercept <- coef_obj$intercept - ret_list$muy +
            drop(ret_list$mux %*% coef_obj$beta)
        coef_obj$beta <- coef_obj$beta * ret_list$scale_x
        return(coef_obj)
    }

    ret_list$unstandardize_coefs <- function(coef_obj) {
        if (!isTRUE(standardize)) {
            return(coef_obj)
        }

        coef_obj$beta <- coef_obj$beta / ret_list$scale_x
        coef_obj$intercept <- coef_obj$intercept + ret_list$muy -
            drop(ret_list$mux %*% coef_obj$beta)
        return(coef_obj)
    }

    return(ret_list)
}

##
## Standardize regression coefficients
##
#' @importFrom Matrix drop
standardize_coefs <- function(intercept, beta, scale_x, mux, muy) {
    return(list(
        intercept = intercept - muy + drop(mux %*% beta),
        beta = beta * scale_x
    ))
}

##
## Unstandardize regression coefficients
##
#' @importFrom Matrix drop
unstandardize_coefs <- function(intercept, beta, scale_x, mux, muy) {
    beta <- beta / scale_x
    return(list(
        intercept = intercept + muy - drop(mux %*% beta),
        beta = beta
    ))
}


## Convenience function for parallel computing
##
## This function handles setting up, using, and closing a potential cluster.
## If no cluster of computing nodes is requested, it will create a proxy to the local
## \code{lapply} function.
##
#' @importFrom parallel clusterEvalQ clusterExport clusterApplyLB stopCluster
#' @importFrom parallel makePSOCKcluster clusterSetRNGStream
#' @importFrom methods is
setupCluster <- function(ncores = 1L, cl = NULL, eval, export, envir = parent.frame()) {
    retlocal <- list(
        lapply = lapply,
        ncores = 1L,
        stopCluster = function() {},
        setSeed = set.seed,
        exportedObject = function(obj) {
            return(obj)
        }
    )

    ret <- retlocal

    ##
    ## Set up a potential cluster
    ##
    if (!is.numeric(ncores) || length(ncores) != 1L || ncores < 1) {
        warning("`ncores` must be a positive integer of length one.")
    } else {
        ret$ncores <- ncores
    }

    tryCatch({
        withCallingHandlers({
            if (is(cl, "cluster") || ncores > 1L) {
                if(!is(cl, "cluster")) {
                    cl <- makePSOCKcluster(ncores)
                    ret$stopCluster <- function() {
                        parallel::stopCluster(cl)
                    }
                }

                ret$ncores <- length(cl)

                if (!missing(eval)) {
                    clusterEvalQ(cl, eval)
                }

                if (!missing(export)) {
                    clusterExport(cl, export, envir = envir)
                }

                ret$setSeed <- function(seed) {
                    clusterSetRNGStream(cl, iseed = seed)
                }

                ret$exportedObject <- function(obj) {
                    substitute(obj)
                }

                ret$lapply <- function(...) {
                    clusterApplyLB(cl, ...)
                }
            }

        }, error = function(...) {
            ret <- retlocal
        })
    }, error = function(e) {
        warning("Error during cluster setup: ", e)
    }, finally = {})

    return(ret)
}

##
## Get the default lambda min ratio, depending on the dimensions of `x`
##
.default_lambda_min_ratio <- function(x) {
    lambda_min_ratio <- tryCatch({
        if (ncol(x) < nrow(x)) {
            1e-4
        } else {
            1e-3
        }
    },
    error = function(...) {
        return(1e-4)
    })
}

##
## Get the EN correction factor
##
#' @useDynLib pense, .registration = TRUE
.en_correction_factor <- function(correction, alpha, lambda) {
    .Call(
        C_en_correction_factor,
        as.integer(correction),
        as.numeric(alpha),
        as.numeric(lambda)
    )
}


