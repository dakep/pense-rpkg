library(testthat)
library(pense)

gradient <- function(coefs, X, y, alpha, lambda) {
    intercept <- coefs[1L]
    beta <- coefs[-1L]

    a <- -crossprod(X, y - X %*% beta - intercept) / length(y)
    p2 <- lambda * (1 - alpha) * beta
    p1 <- lambda * alpha * sign(beta)

    gr <- drop(a + p1 + p2)

    ## If the i-th element of beta is 0, the subgradient of
    ## the L1 norm at the i-th element is [-1; 1]
    gr[(abs(gr) <= 1) & (abs(beta) < .Machine$double.eps)] <- 0

    return(gr)
}

en_options <- en_options_dal(eps = 1e-9);
EQUALITY_TOLERANCE = 1e-5
PROTECTED_VARS <- c("gradient", "en_options", "augment",
                    "EQUALITY_TOLERANCE", "PROTECTED_VARS")

test_that("LASSO", {
    ##
    ## A fairly simple case
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda1 <- 0.02

    enres <- elnet(X, y, 1, lambda1, intercept = FALSE, options = en_options)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, 1, lambda1),
                 numeric(p))

    ## check if result is the same as for lars
    if (requireNamespace("lars", quietly = TRUE)) {
        larsobj <- lars::lars(X, y, type = "lasso", normalize = FALSE,
                              intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = n * lambda1, mode = "lambda")

        expect_equal(enres$coefficients, c(0, larsres),
                     tolerance = EQUALITY_TOLERANCE)
    }

    remove(list = setdiff(ls(), PROTECTED_VARS))

    ##
    ## A fairly simple case with large X values
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- 100 * matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda1 <- 0.02

    enres <- elnet(X, y, 1, lambda1, intercept = FALSE, options = en_options)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, 1, lambda1),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    ## check if result is the same as for lars
    if (requireNamespace("lars", quietly = TRUE)) {
        larsobj <- lars::lars(X, y, type = "lasso", normalize = FALSE, intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = n * lambda1, mode = "lambda")

        expect_equal(enres$coefficients, c(0, larsres),
                     tolerance = EQUALITY_TOLERANCE)
    }

    remove(list = setdiff(ls(), PROTECTED_VARS))

    ##
    ## Some more observations
    ##
    n <- 2000L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda1 <- 0.2

    enres <- elnet(X, y, 1, lambda1, intercept = FALSE, options = en_options)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, 1, lambda1),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    ## check if result is the same as for lars
    if (requireNamespace("lars", quietly = TRUE)) {
        larsobj <- lars::lars(X, y, type = "lasso", normalize = FALSE,
                              intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = n * lambda1, mode = "lambda")

        expect_equal(enres$coefficients, c(0, larsres),
                     tolerance = EQUALITY_TOLERANCE)
    }

    remove(list = setdiff(ls(), PROTECTED_VARS))

    ##
    ## More observations than variables with reasonable regularization
    ##
    n <- 100L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda1 <- 0.02

    enres <- elnet(X, y, 1, lambda1, intercept = FALSE, options = en_options)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, 1, lambda1),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    ## check if result is the same as for lars
    if (requireNamespace("lars", quietly = TRUE)) {
        larsobj <- lars::lars(X, y, type = "lasso", normalize = FALSE,
                              intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = n * lambda1, mode = "lambda")

        expect_equal(enres$coefficients, c(0, larsres),
                     tolerance = EQUALITY_TOLERANCE)
    }

    remove(list = setdiff(ls(), PROTECTED_VARS))

    ##
    ## More observations than variables with almost no regularization
    ##
    n <- 100L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda1 <- 0.002

    enres <- elnet(X, y, 1, lambda1, intercept = FALSE, options = en_options)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, 1, lambda1),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    ## check if result is the same as for lars
    if (requireNamespace("lars", quietly = TRUE)) {
        larsobj <- lars::lars(X, y, type = "lasso", normalize = FALSE,
                              intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = n * lambda1, mode = "lambda")

        expect_equal(enres$coefficients, c(0, larsres),
                     tolerance = EQUALITY_TOLERANCE)
    }

    remove(list = setdiff(ls(), PROTECTED_VARS))
})

test_that("Ridge", {
    augment <- function(X, y, lambda2, leading1s = TRUE) {
        d <- dim(X)

        ext <- diag(sqrt(n * lambda2), ncol = d[2L], nrow = d[2L])

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

    enres <- elnet(X, y, 0, lambda2, intercept = TRUE, options = en_options)
    olsres <- .lm.fit(au$X, au$y)

    expect_equal(enres$coefficients, olsres$coefficients,
                 tolerance = EQUALITY_TOLERANCE)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, 0, lambda2),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    remove(list = setdiff(ls(), PROTECTED_VARS))

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

    enres <- elnet(X, y, 0, lambda2, intercept = TRUE, options = en_options)
    olsres <- .lm.fit(au$X, au$y)

    expect_equal(enres$coefficients, olsres$coefficients,
                 tolerance = EQUALITY_TOLERANCE)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, 0, lambda2),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    remove(list = setdiff(ls(), PROTECTED_VARS))

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

    enres <- elnet(X, y, 0, lambda2, intercept = TRUE, options = en_options)
    olsres <- .lm.fit(au$X, au$y)

    expect_equal(enres$coefficients, olsres$coefficients,
                 tolerance = EQUALITY_TOLERANCE)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, 0, lambda2),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    remove(list = setdiff(ls(), PROTECTED_VARS))

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

    enres <- elnet(X, y, 0, lambda2, intercept = TRUE, options = en_options)
    olsres <- .lm.fit(au$X, au$y)

    expect_equal(enres$coefficients, olsres$coefficients,
                 tolerance = EQUALITY_TOLERANCE)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, 0, lambda2),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    remove(list = setdiff(ls(), PROTECTED_VARS))
})

test_that("EN", {
    augment <- function(X, y, lambda2, leading1s = TRUE) {
        d <- dim(X)

        ext <- diag(sqrt(lambda2 * n * 2), ncol = d[2L], nrow = d[2L])

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
    ## A fairly simple case without 2 norm penalty
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda2 <- 0
    lambda1 <- 0.05

    lambda <- 2 * lambda2 + lambda1
    alpha <- lambda1 / (2 * lambda2 + lambda1)

    enres <- elnet(X, y, alpha, lambda, intercept = FALSE, options = en_options)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, alpha, lambda),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    ## check if we can match the result by augmenting the data
    au <- augment(X, y, lambda2, leading1s = FALSE)

    elau <- elnet(au$X, au$y, alpha = 1, n * lambda1 / (n + p),
                  intercept = FALSE, options = en_options)
    expect_equal(enres$coefficients[-1L], elau$coefficients[-1L],
                 tolerance = EQUALITY_TOLERANCE)

    ## check if results match with lars
    if (requireNamespace("lars", quietly = TRUE)) {
        larsobj <- lars::lars(au$X, au$y, type = "lasso", normalize = FALSE,
                              intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = n * lambda1, mode = "lambda")
        expect_equal(enres$coefficients[-1L], larsres,
                     tolerance = EQUALITY_TOLERANCE)
    }

    remove(list = setdiff(ls(), PROTECTED_VARS))

    ##
    ## A fairly simple case
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda2 <- 0.02
    lambda1 <- 0.1

    lambda <- 2 * lambda2 + lambda1
    alpha <- lambda1 / (2 * lambda2 + lambda1)

    enres <- elnet(X, y, alpha, lambda, intercept = FALSE, options = en_options)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, alpha, lambda),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    ## check if we can match the result by augmenting the data
    au <- augment(X, y, lambda2, leading1s = FALSE)

    elau <- elnet(au$X, au$y, alpha = 1, n * lambda1 / (n + p),
                  intercept = FALSE, options = en_options)
    expect_equal(enres$coefficients[-1L], elau$coefficients[-1L],
                 tolerance = EQUALITY_TOLERANCE)

    ## check if results match with lars
    if (requireNamespace("lars", quietly = TRUE)) {
        larsobj <- lars::lars(au$X, au$y, type = "lasso", normalize = FALSE,
                              intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = n * lambda1, mode = "lambda")
        expect_equal(enres$coefficients[-1L], larsres,
                     tolerance = EQUALITY_TOLERANCE)
    }

    remove(list = setdiff(ls(), PROTECTED_VARS))

    ##
    ## A fairly simple case with many observations
    ##
    n <- 2000L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    lambda2 <- 0.005
    lambda1 <- 0.02

    lambda <- 2 * lambda2 + lambda1
    alpha <- lambda1 / (2 * lambda2 + lambda1)

    enres <- elnet(X, y, alpha, lambda, intercept = FALSE, options = en_options)

    ## check if zero is in the gradient
    expect_equal(gradient(enres$coefficients, X, y, alpha, lambda),
                 numeric(p), tolerance = EQUALITY_TOLERANCE)

    ## check if we can match the result by augmenting the data
    au <- augment(X, y, lambda2, leading1s = FALSE)

    elau <- elnet(au$X, au$y, alpha = 1, n * lambda1 / (n + p),
                  intercept = FALSE, options = en_options)
    expect_equal(enres$coefficients[-1L], elau$coefficients[-1L],
                 tolerance = EQUALITY_TOLERANCE)

    ## check if results match with lars
    if (requireNamespace("lars", quietly = TRUE)) {
        larsobj <- lars::lars(au$X, au$y, type = "lasso", normalize = FALSE,
                              intercept = FALSE)
        larsres <- lars::coef.lars(larsobj, s = n * lambda1, mode = "lambda")
        expect_equal(enres$coefficients[-1L], larsres,
                     tolerance = EQUALITY_TOLERANCE)
    }

    remove(list = setdiff(ls(), PROTECTED_VARS))
})

test_that("EN - Bugs", {
    n <- 100L
    p <- 100L
    ##
    ## Test behaviour when no observations are given
    ##
    X <- matrix(0.0, ncol = p, nrow = 0L)
    y <- numeric(0L)
    res <- elnet(X, y, alpha = 0.5, lambda = 2, options = en_options)

    # All coefficients should be zero and residuals of length zero
    expect_identical(res$coefficients, numeric(p + 1L))
    expect_identical(res$residuals, numeric(0L))

    ##
    ## Test behaviour when no columns are given
    ##
    X <- matrix(0.0, ncol = 0L, nrow = n)
    y <- rnorm(n)
    res <- elnet(X, y, alpha = 0.5, lambda = 2, addLeading1s = FALSE,
                 options = en_options)

    # All coefficients should be zero and residuals of length zero
    expect_identical(res$coefficients, numeric(0L))
    expect_identical(res$residuals, y)

    ##
    ## Test behaviour when only the column of 1's is given (i.e., average)
    ##
    X <- matrix(0.0, ncol = 0L, nrow = n)
    y <- rnorm(n)
    res <- elnet(X, y, alpha = 0.5, lambda = 2, options = en_options)

    # All coefficients should be zero and residuals of length zero
    expect_equal(res$coefficients, mean(y))
    expect_equal(res$residuals, y - mean(y))
})
