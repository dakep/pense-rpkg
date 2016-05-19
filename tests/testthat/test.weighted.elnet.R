.elnet.wfit.rimpl <- function(X, y, weights, alpha, lambda) {
    k <- sqrt(weights)
    Xweight <- diag(k) %*% X
    yweight <- k * y

    etat <- crossprod(k, Xweight) / sum(weights)
    Xweightorth <- Xweight - k %*% etat

    res <- elnet(Xweightorth, yweight, alpha = alpha, lambda = lambda, centering = FALSE)
    res$residuals <- drop(y - X %*% res$coefficients[-1L])
    res$coefficients[1L] <- weighted.mean(res$residuals, weights)
    res$residuals <- res$residuals - res$coefficients[1L]

    return(res)
}

test_that("Weighted LASSO", {
    alpha <- 1
    lambda <- 0.02

    ##
    ## A fairly simple case
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)

    ##
    ## A fairly simple case with large X values
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- 100 * matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)

    ##
    ## Some more observations
    ##
    n <- 2000L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)

    ##
    ## More observations than variables with reasonable regularization
    ##
    n <- 100L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)

    ##
    ## More observations than variables with almost no regularization
    ##
    n <- 100L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)
})

test_that("Weighted Ridge", {
    alpha <- 0
    lambda <- 0.02

    ##
    ## A fairly simple case
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)

    ##
    ## Some more observations
    ##
    n <- 2000L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)

    ##
    ## More observations than variables with reasonable regularization
    ##
    n <- 100L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)

    ##
    ## More observations than variables with almost no regularization
    ##
    n <- 100L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)

})

test_that("Weighted EN", {
    alpha <- 0.5
    lambda <- 0.01

    ##
    ## A fairly simple case without 2 norm penalty
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)


    ##
    ## A fairly simple case
    ##
    n <- 200L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)


    ##
    ## A fairly simple case with many observations
    ##
    n <- 2000L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)

    ##
    ## More observations than variables with almost no regularization
    ##
    n <- 100L
    p <- 150L

    set.seed(1234)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    weights <- rchisq(n, 1)

    enres <- elnet(X, y, weights = weights, alpha = alpha, lambda = lambda)
    enres.rimpl <- .elnet.wfit.rimpl(X, y, weights, alpha, lambda)

    expect_equivalent(enres.rimpl$coefficients, enres$coefficients)
})
