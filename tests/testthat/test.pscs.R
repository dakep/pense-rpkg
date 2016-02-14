test_that("Compute PSCs", {
    pscs_R <- function(X, y, intercept = TRUE) {
        if (identical(intercept, TRUE)) {
            X <- cbind(1, X)
        }

        XextCp <- crossprod(X)
        Xsqrt <- chol(XextCp)
        tXsqrtXX <- backsolve(Xsqrt, t(X), transpose = TRUE)

        H <- X %*% backsolve(Xsqrt, tXsqrtXX)
        resids <- as.vector(y - H %*% y)
        Wd <- resids / (1 - diag(H))

        G <- tXsqrtXX %*% diag(Wd)
        Q <- tcrossprod(G)

        egQ <- eigen(Q)
        # Z <- apply(egQ$vectors, 2, crossprod, tXsqrtXX)
        Z <- crossprod(tXsqrtXX, egQ$vectors)
        return(Z[ , order(egQ$values)])
    }

    ##
    ## Normal use case
    ##
    n <- 200L
    p <- 50L
    set.seed(12345)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    res_r <- pscs_R(X, y)
    res_c <- penseinit::prinsens(X, y)

    expect_lt(mean(abs(abs(res_r / res_c) - 1)), 1e-9)

    res_r <- pscs_R(X, y, intercept = FALSE)
    res_c <- penseinit::prinsens(X, y, intercept = FALSE)

    expect_lt(mean(abs(abs(res_r / res_c) - 1)), 1e-9)

    ##
    ## Singular matrix
    ##
    n <- 20L
    p <- 50L
    set.seed(12345)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    expect_error(penseinit::prinsens(X, y), regexp = "singular")
})
