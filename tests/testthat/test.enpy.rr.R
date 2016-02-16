test_that("enpy.rr", {
    source("enpy.old.R", local = TRUE)

    n <- 500L
    p <- 50L
    set.seed(12345)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    ##
    ## No penalization at all
    ##
    lambda1 <- 0
    lambda2 <- 0

    new <- penseinit::enpy(X, y, lambda1 = lambda1, lambda2 = lambda2,
                           deltaesc = 0.5, cc.scale = 1.54764, psc.method = "rr",
                           prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
                           py.nit = 5, en.tol = 1e-8)

    target <- enpy(X, y, lambda1 = lambda1, lambda2 = lambda2,
                   deltaesc = 0.5, cc.scale = 1.54764, psc.method = "rr",
                   prosac = 0.8, clean.method = "proportion", prop = 0.4,
                   py.nit = 5, en.tol = 1e-8)

    diffs <- apply(new$coeff, 2, function(x) {
        min(colSums(abs(target$coeff - x)))
    })

    expect_lt(sum(abs(diffs)), 1e-9)


    ##
    ## Ridge regression penalty
    ##
    lambda1 <- 0
    lambda2 <- 1.5

    new <- penseinit::enpy(X, y, lambda1 = lambda1, lambda2 = lambda2,
                           deltaesc = 0.5, cc.scale = 1.54764, psc.method = "rr",
                           prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
                           py.nit = 5, en.tol = 1e-8)

    target <- enpy(X, y, lambda1 = lambda1, lambda2 = lambda2,
                   deltaesc = 0.5, cc.scale = 1.54764, psc.method = "rr",
                   prosac = 0.8, clean.method = "proportion", prop = 0.4,
                   py.nit = 5, en.tol = 1e-8)

    diffs <- apply(new$coeff, 2, function(x) {
        min(colSums(abs(target$coeff - x)))
    })

    expect_lt(sum(abs(diffs)), 1e-9)
})
