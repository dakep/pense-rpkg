## Generate a grid of lambda values
## The grid is *INDEPENDENT* from the sample size
##
#' @importFrom stats mad median
#' @importFrom robustbase scaleTau2 covGK Mwgt
.build_lambda_grid_s <- function(x, y, alpha, nlambda, lambda_min_ratio, options) {
    # Compute intercept and scale estimate in intercept-only model
    loc_scale_est <- .mlocscale_simultaneous(
        y,
        delta = options$bdp,
        c1 = options$cc,
        c2 = options$cc
    )

    resid <- (y - loc_scale_est["location"])
    ## Old code pense <= 1.0.8
    # w_0 <- Mwgt(resid, cc = options$cc * loc_scale_est["scale"], "bisquare")
    # ws_0 <- w_0 / mean(w_0)
    ## New code pense > 1.0.8
    drho_wgt <- Mwgt(
        resid,
        cc = loc_scale_est["scale"] * options$cc,
        psi = "bisquare"
    )
    ws_0 <- drho_wgt * (2 * loc_scale_est["scale"]^2 / mean(drho_wgt * resid^2))
    if (is.null(dim(x))) {
        abs_grad <- abs(mean(x * resid * ws_0))
    }
    abs_grad <- abs(colMeans(x * resid * ws_0))

    first_var <- which.max(abs_grad)
    first_var_grad <- abs_grad[[first_var]]

    lmax <- first_var_grad / max(0.01, alpha)
    lmin <- lambda_min_ratio * lmax

    d <- diff(log(c(lmin, lmax))) / (nlambda - 2L)
    lambda <- exp(seq(
        log(lmin),
        log(lmax) + d,
        by = d)
    )

    return(lambda)
}

#' @importFrom robustbase Mwgt MrhoInf
.lambda_max_m <- function(x, y, alpha, scale_init, bdp, cc) {
    y <- as.numeric(y)
    # Compute intercept and scale estimate in intercept-only model
    loc_scale_est <- .mloc(
        y,
        delta = bdp,
        c1 = cc,
        scale_init = scale_init
    )

    resid <- (y - loc_scale_est["location"])
    ## Old code pense <= 1.0.8
    # w_0 <- Mwgt(resid, cc = options$cc * scale_init, "bisquare")
    # ws_0 <- w_0 / mean(w_0)
    ## New code pense > 1.0.8
    ws_0 <- Mwgt(
        resid,
        cc = scale_init * cc,
        psi = "bisquare"
    ) / MrhoInf(scale_init * cc, "bisquare")

    abs_grad <- if (is.null(dim(x))) {
        abs(mean(x * resid * ws_0))
    } else {
        abs(colMeans(x * resid * ws_0))
    }

    first_var <- which.max(abs_grad)
    first_var_grad <- abs_grad[[first_var]]

    return(first_var_grad / max(0.001, alpha))
}

#' @importFrom robustbase Mwgt Mchi
#' @importFrom stats weighted.mean
.mlocscale_simultaneous <- function (x, c1, delta, c2, eps = 1e-8, maxit = 200) {
    loc <- median(x)
    scale <- mad(x)

    conv_tol <- eps * scale
    prev_loc <- 0
    prev_scale <- 1
    it <- 0L
    repeat {
        r <- x - loc
        w_loc <- Mwgt(r, c1 * scale, "bisquare")
        w_scale <- Mchi(r, scale * c2, "bisquare")
        prev_loc <- loc
        prev_scale <- scale

        loc <- weighted.mean(x, w_loc)
        scale <- prev_scale * sqrt(mean(w_scale)) / sqrt(delta)

        if ((abs(prev_loc - loc) < conv_tol) && (abs(prev_scale - scale) < conv_tol)) {
            break
        }

        it <- it + 1L
        if (it >= maxit) {
            break
        }
    }

    structure(
        c(location = loc, scale = scale),
        it = it
    )
}

#' @importFrom robustbase Mwgt
#' @importFrom stats weighted.mean
.mloc <- function (x, c1, delta, eps = 1e-8, maxit = 200, scale_init) {
    loc <- median(x)

    c1 <- scale_init * c1
    conv_tol <- eps * scale_init
    prev_loc <- 0
    it <- 0L
    repeat {
        w <- Mwgt(x - loc, c1, "bisquare")
        prev_loc <- loc
        loc <- weighted.mean(x, w)
        if (is.na(loc)) {
            # All weights are 0
            loc <- prev_loc
        }
        if (abs(prev_loc - loc) < conv_tol) {
            break
        }

        it <- it + 1L
        if (it >= maxit) {
            break
        }
    }

    structure(
        c(location = loc, scale = mscale(as.numeric(x - loc), delta = delta)),
        it = it
    )
}

