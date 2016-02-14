test_that("LASSO", {
    if (requireNamespace("lars", quietly = TRUE)) {
        ##
        ## A fairly simple case
        ##
        n <- 200L
        p <- 150L

        set.seed(1234)
        X <- matrix(rnorm(n * p), ncol = p)
        y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

        lambda1 <- 0.2

        enres <- elnet(X, y, 1, lambda1 / n, eps = 1e-10, centering = FALSE)

        larsobj <- lars::lars(X, y, type = "lasso", normalize = FALSE, intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = lambda1, mode = "lambda")

        expect_equal(enres$coefficients, c(0, larsres))

        remove(list = ls())

        ##
        ## A fairly simple case with large X values
        ##
        n <- 200L
        p <- 150L

        set.seed(1234)
        X <- 100 * matrix(rnorm(n * p), ncol = p)
        y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

        lambda1 <- 0.2

        enres <- elnet(X, y, 1, lambda1 / n, eps = 1e-10, centering = FALSE)

        larsobj <- lars::lars(X, y, type = "lasso", normalize = FALSE, intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = lambda1, mode = "lambda")

        expect_equal(enres$coefficients, c(0, larsres))

        remove(list = ls())

        ##
        ## Some more observations
        ##
        n <- 2000L
        p <- 150L

        set.seed(1234)
        X <- matrix(rnorm(n * p), ncol = p)
        y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

        lambda1 <- 0.2

        enres <- elnet(X, y, 1, lambda1 / n, eps = 1e-10, centering = FALSE)

        larsobj <- lars::lars(X, y, type = "lasso", normalize = FALSE, intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = lambda1, mode = "lambda")

        expect_equal(enres$coefficients, c(0, larsres))

        remove(list = ls())

        ##
        ## More observations than variables with reasonable regularization
        ##
        n <- 100L
        p <- 150L

        set.seed(1234)
        X <- matrix(rnorm(n * p), ncol = p)
        y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

        lambda1 <- 3

        enres <- elnet(X, y, 1, lambda1 / n, eps = 1e-10, centering = FALSE)

        larsobj <- lars::lars(X, y, type = "lasso", normalize = FALSE, intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = lambda1, mode = "lambda")

        expect_equal(enres$coefficients, c(0, larsres))

        remove(list = ls())

        ##
        ## More observations than variables with almost no regularization
        ##
        n <- 100L
        p <- 150L

        set.seed(1234)
        X <- matrix(rnorm(n * p), ncol = p)
        y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

        lambda1 <- 1

        enres <- elnet(X, y, 1, lambda1 / n, eps = 1e-10, maxit = 1e5, centering = FALSE)

        larsobj <- lars::lars(X, y, type = "lasso", normalize = FALSE, intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = lambda1, mode = "lambda")

        expect_equal(enres$coefficients, c(0, larsres), tolerance = .Machine$double.eps ^ 0.3)

        remove(list = ls())
    }
})

test_that("Ridge", {
    augment <- function(X, y, lambda2, leading1s = TRUE) {
        d <- dim(X)

        ext <- diag(sqrt(lambda2), ncol = d[2L], nrow = d[2L])

        if (identical(leading1s, TRUE)) {
            X <- cbind(1, X)
            ext <- cbind(0, ext)
        }

        return(list(
            X = rbind(X, ext),
            y = c(y, numeric(d[2L]))
        ))
    }

    ##
    ## A fairly simple case
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda2 <- 0.2

    au <- augment(X, y, lambda2)

    enres <- elnet(X, y, 0, lambda2 / n, eps = 1e-10, centering = TRUE)
    olsres <- .lm.fit(au$X, au$y)

    expect_equal(enres$coefficients, olsres$coefficients)

    ##
    ## Some more observations
    ##
    n <- 2000L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda2 <- 0.2
    au <- augment(X, y, lambda2)

    enres <- elnet(X, y, 0, lambda2 / n, eps = 1e-10, centering = TRUE)
    olsres <- .lm.fit(au$X, au$y)

    expect_equal(enres$coefficients, olsres$coefficients)

    ##
    ## More observations than variables with reasonable regularization
    ##
    n <- 100L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda2 <- 3
    au <- augment(X, y, lambda2)

    enres <- elnet(X, y, 0, lambda2 / n, eps = 1e-10, centering = TRUE)
    olsres <- .lm.fit(au$X, au$y)

    expect_equal(enres$coefficients, olsres$coefficients)

    ##
    ## More observations than variables with almost no regularization
    ##
    n <- 100L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda2 <- 1
    au <- augment(X, y, lambda2)

    enres <- elnet(X, y, 0, lambda2 / n, eps = 1e-10, maxit = 1e5, centering = TRUE)
    olsres <- .lm.fit(au$X, au$y)

    expect_equal(enres$coefficients, olsres$coefficients)
})

test_that("EN", {
    augment <- function(X, y, lambda2, leading1s = TRUE) {
        d <- dim(X)

        ext <- diag(sqrt(lambda2), ncol = d[2L], nrow = d[2L])

        if (identical(leading1s, TRUE)) {
            X <- cbind(1, X)
            ext <- cbind(0, ext)
        }

        return(list(
            X = rbind(X, ext),
            y = c(y, numeric(d[2L]))
        ))
    }

    ##
    ## A fairly simple case
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda2 <- 0.5
    lambda1 <- 1

    lambda <- 2 * lambda2 + lambda1
    alpha <- lambda1 / (2 * lambda2 + lambda1)

    enres <- elnet(X, y, alpha, lambda / n, eps = 1e-10, centering = FALSE)

    au <- augment(X, y, 2 * lambda2, leading1s = FALSE)

    larsobj <- lars::lars(au$X, au$y, type = "lasso", normalize = FALSE, intercept = FALSE)
    larsres <- lars::coef.lars(larsobj, s = lambda1, mode = "lambda")

    elau <- elnet(au$X, au$y, alpha = 1, lambda / (2 * (n + p)), eps = 1e-10, centering = FALSE)

    expect_equal(enres$coefficients[-1L], larsres)
    expect_equal(enres$coefficients[-1L], elau$coefficients[-1L])
})
