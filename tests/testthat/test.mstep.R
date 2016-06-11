test_that("mstep", {
    ##
    ## Test case 1
    ##
    n <- 500L
    p <- 50L
    set.seed(12345)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)


    lambda <- 0.01
    cc <- 1.5
    init.scale <- 1.1
    init.coef <- numeric(p + 1)
    ctrl <- pense.control()

    ## LASSO penalty
    alpha <- 1

    mr.r <- pense:::pensemstep.rimpl(X, y, cc = cc,
                                         init.scale = init.scale, init.coef = init.coef,
                                         alpha = alpha, lambda = lambda, control = ctrl)

    mr.c <- pense:::pensemstep(X, y, cc = cc,
                                   init.scale = init.scale, init.coef = init.coef,
                                   alpha = alpha, lambda = lambda, control = ctrl)

    expect_equal(mr.c$intercept, mr.r$intercept)
    expect_equal(mr.c$beta, mr.r$beta)
    expect_equal(mr.c$objF, mr.r$objF)


    ## EN penalty
    alpha <- 0.5

    mr.r <- pense:::pensemstep.rimpl(X, y, cc = cc,
                                         init.scale = init.scale, init.coef = init.coef,
                                         alpha = alpha, lambda = lambda, control = ctrl)

    mr.c <- pense:::pensemstep(X, y, cc = cc,
                                   init.scale = init.scale, init.coef = init.coef,
                                   alpha = alpha, lambda = lambda, control = ctrl)

    expect_equal(mr.c$intercept, mr.r$intercept)
    expect_equal(mr.c$beta, mr.r$beta)
    expect_equal(mr.c$objF, mr.r$objF)

    ##
    ## Test case 2
    ##
    remove(list = ls())
    n <- 500L
    p <- 50L
    set.seed(12345)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)


    lambda <- 0.01
    cc <- 1.5
    init.scale <- 1.1
    init.coef <- rep.int(10, p + 1L)
    ctrl <- pense.control()

    ## LASSO penalty
    alpha <- 1

    mr.r <- pense:::pensemstep.rimpl(X, y, cc = cc,
                                         init.scale = init.scale, init.coef = init.coef,
                                         alpha = alpha, lambda = lambda, control = ctrl)

    mr.c <- pense:::pensemstep(X, y, cc = cc,
                                   init.scale = init.scale, init.coef = init.coef,
                                   alpha = alpha, lambda = lambda, control = ctrl)

    expect_equal(mr.c$intercept, mr.r$intercept)
    expect_equal(mr.c$beta, mr.r$beta)
    expect_equal(mr.c$objF, mr.r$objF)


    ## EN penalty
    alpha <- 0.5

    mr.r <- pense:::pensemstep.rimpl(X, y, cc = cc,
                                         init.scale = init.scale, init.coef = init.coef,
                                         alpha = alpha, lambda = lambda, control = ctrl)

    mr.c <- pense:::pensemstep(X, y, cc = cc,
                                   init.scale = init.scale, init.coef = init.coef,
                                   alpha = alpha, lambda = lambda, control = ctrl)

    expect_equal(mr.c$intercept, mr.r$intercept)
    expect_equal(mr.c$beta, mr.r$beta)
    expect_equal(mr.c$objF, mr.r$objF)

    ##
    ## Test case 3
    ##
    remove(list = ls())
    n <- 500L
    p <- 100L
    set.seed(12345)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)


    lambda <- 0.02
    cc <- 1.5
    init.scale <- 1.1
    init.coef <- rep.int(4, p + 1L)
    ctrl <- pense.control()

    ## LASSO penalty
    alpha <- 1

    mr.r <- pense:::pensemstep.rimpl(X, y, cc = cc,
                                         init.scale = init.scale, init.coef = init.coef,
                                         alpha = alpha, lambda = lambda, control = ctrl)

    mr.c <- pense:::pensemstep(X, y, cc = cc,
                                   init.scale = init.scale, init.coef = init.coef,
                                   alpha = alpha, lambda = lambda, control = ctrl)

    expect_equal(mr.c$intercept, mr.r$intercept)
    expect_equal(mr.c$beta, mr.r$beta)
    expect_equal(mr.c$objF, mr.r$objF)


    ## EN penalty
    alpha <- 0.5

    mr.r <- pense:::pensemstep.rimpl(X, y, cc = cc,
                                         init.scale = init.scale, init.coef = init.coef,
                                         alpha = alpha, lambda = lambda, control = ctrl)

    mr.c <- pense:::pensemstep(X, y, cc = cc,
                                   init.scale = init.scale, init.coef = init.coef,
                                   alpha = alpha, lambda = lambda, control = ctrl)

    expect_equal(mr.c$intercept, mr.r$intercept)
    expect_equal(mr.c$beta, mr.r$beta)
    expect_equal(mr.c$objF, mr.r$objF)

    ##
    ## Test case 4 -- should not converge
    ##
    remove(list = ls())
    n <- 500L
    p <- 100L
    set.seed(12345)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)


    lambda <- 0.02
    cc <- 1.5
    init.scale <- 1.1
    init.coef <- rep.int(4, p + 1L)
    ctrl <- pense.control(pense.maxit = 2)

    ## LASSO penalty
    alpha <- 1

    expect_warning(
        pense:::pensemstep.rimpl(X, y, cc = cc,
                                     init.scale = init.scale, init.coef = init.coef,
                                     alpha = alpha, lambda = lambda, control = ctrl),
        regexp = 'did not converge'
    )

    expect_warning(
        pense:::pensemstep(X, y, cc = cc,
                               init.scale = init.scale, init.coef = init.coef,
                               alpha = alpha, lambda = lambda, control = ctrl),
        regexp = 'did not converge'
    )


    ##
    ## Test case 5 -- leverage points
    ##
    remove(list = ls())
    n <- 500L
    p <- 50L
    set.seed(12345)
    X <- matrix(rnorm(n * p), ncol = p)
    y <- 2 + X %*% c(1, 1, 1, rep.int(0, p - 3L)) + rnorm(n)

    X[1:50, ] <- 200 + X[1:50, ]

    lambda <- 0.01
    cc <- 1.5
    init.scale <- 1.1
    init.coef <- numeric(p + 1)
    ctrl <- pense.control()

    ## LASSO penalty
    alpha <- 1

    mr.r <- pense:::pensemstep.rimpl(X, y, cc = cc,
                                         init.scale = init.scale, init.coef = init.coef,
                                         alpha = alpha, lambda = lambda, control = ctrl)

    mr.c <- pense:::pensemstep(X, y, cc = cc,
                                   init.scale = init.scale, init.coef = init.coef,
                                   alpha = alpha, lambda = lambda, control = ctrl)

    expect_equal(mr.c$intercept, mr.r$intercept)
    expect_equal(mr.c$beta, mr.r$beta)
    expect_equal(mr.c$objF, mr.r$objF)


    ## EN penalty
    alpha <- 0.5

    mr.r <- pense:::pensemstep.rimpl(X, y, cc = cc,
                                         init.scale = init.scale, init.coef = init.coef,
                                         alpha = alpha, lambda = lambda, control = ctrl)

    mr.c <- pense:::pensemstep(X, y, cc = cc,
                                   init.scale = init.scale, init.coef = init.coef,
                                   alpha = alpha, lambda = lambda, control = ctrl)

    expect_equal(mr.c$intercept, mr.r$intercept)
    expect_equal(mr.c$beta, mr.r$beta)
    expect_equal(mr.c$objF, mr.r$objF)

})
