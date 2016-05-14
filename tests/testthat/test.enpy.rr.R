test_that("enpy.rr", {
    if (requireNamespace("lars", quietly = TRUE) &&
        requireNamespace("robustbase", quietly = TRUE)) {
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

        expect_lt(sum(abs(diffs)), 1e-8)


        ##
        ## Ridge regression penalty
        ##
        lambda <- 0.1

        # alpha <- 0
        # rr.rr.old <- penseinit::enpy(X, y, lambda1 = alpha * lambda, lambda2 = 0.5 * (1 - alpha) * lambda,
        #                           deltaesc = 0.5, cc.scale = 1.54764, psc.method = "rr",
        #                           prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
        #                           py.nit = 5, en.tol = 1e-8)
        #
        # rr.mn.old <- penseinit::enpy(X, y, lambda1 = alpha * lambda, lambda2 = 0.5 * (1 - alpha) * lambda,
        #                           deltaesc = 0.5, cc.scale = 1.54764, psc.method = "Mn",
        #                           prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
        #                           py.nit = 5, en.tol = 1e-8)
        #
        # alpha <- 1
        # lr.rr.old <- penseinit::enpy(X, y, lambda1 = alpha * lambda, lambda2 = 0.5 * (1 - alpha) * lambda,
        #                           deltaesc = 0.5, cc.scale = 1.54764, psc.method = "rr",
        #                           prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
        #                           py.nit = 5, en.tol = 1e-8)
        # lr.mn.old <- penseinit::enpy(X, y, lambda1 = alpha * lambda, lambda2 = 0.5 * (1 - alpha) * lambda,
        #                           deltaesc = 0.5, cc.scale = 1.54764, psc.method = "Mn",
        #                           prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
        #                           py.nit = 5, en.tol = 1e-8)
        #
        # alpha <- 0.5
        # en.rr.old <- penseinit::enpy(X, y, lambda1 = alpha * lambda, lambda2 = 0.5 * (1 - alpha) * lambda,
        #                           deltaesc = 0.5, cc.scale = 1.54764, psc.method = "rr",
        #                           prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
        #                           py.nit = 5, en.tol = 1e-8)
        # en.mn.old <- penseinit::enpy(X, y, lambda1 = alpha * lambda, lambda2 = 0.5 * (1 - alpha) * lambda,
        #                           deltaesc = 0.5, cc.scale = 1.54764, psc.method = "Mn",
        #                           prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
        #                           py.nit = 5, en.tol = 1e-8)
        #
        #
        # save(rr.rr.old, rr.mn.old,
        #      lr.rr.old, lr.mn.old,
        #      en.rr.old, en.mn.old,
        #      file = "~/Desktop/tmp.rda")

        lambda <- 0.1

        alpha <- 0
        rr.rr.new <- penseinit::enpy(X, y, alpha = alpha, lambda = lambda,
                                     deltaesc = 0.5, cc.scale = 1.54764, psc.method = "rr",
                                     prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
                                     py.nit = 5, en.tol = 1e-8)

        rr.mn.new <- penseinit::enpy(X, y, alpha = alpha, lambda = lambda,
                                     deltaesc = 0.5, cc.scale = 1.54764, psc.method = "Mn",
                                     prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
                                     py.nit = 5, en.tol = 1e-8)

        alpha <- 1
        lr.rr.new <- penseinit::enpy(X, y, alpha = alpha, lambda = lambda,
                                     deltaesc = 0.5, cc.scale = 1.54764, psc.method = "rr",
                                     prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
                                     py.nit = 5, en.tol = 1e-8)
        lr.mn.new <- penseinit::enpy(X, y, alpha = alpha, lambda = lambda,
                                     deltaesc = 0.5, cc.scale = 1.54764, psc.method = "Mn",
                                     prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
                                     py.nit = 5, en.tol = 1e-8)

        alpha <- 0.5
        en.rr.new <- penseinit::enpy(X, y, alpha = alpha, lambda = lambda,
                                     deltaesc = 0.5, cc.scale = 1.54764, psc.method = "rr",
                                     prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
                                     py.nit = 5, en.tol = 1e-8)
        en.mn.new <- penseinit::enpy(X, y, alpha = alpha, lambda = lambda,
                                     deltaesc = 0.5, cc.scale = 1.54764, psc.method = "Mn",
                                     prosac = 0.8, clean.method = "proportion", prop = 0.4, C.res = NULL,
                                     py.nit = 5, en.tol = 1e-8)


        load("~/Desktop/tmp.rda")


        sum(abs(rr.rr.new$coeff - rr.rr.old$coeff))
        sum(abs(rr.mn.new$coeff - rr.mn.old$coeff))

        sum(abs(lr.rr.new$coeff - lr.rr.old$coeff))
        sum(abs(lr.mn.new$coeff - lr.mn.old$coeff))

        sum(abs(en.rr.new$coeff - en.rr.old$coeff))
        sum(abs(en.mn.new$coeff - en.mn.old$coeff))
    }
})
