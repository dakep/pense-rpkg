test_that("pen.s.reg", {
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

    pr.r <- pense:::pen.s.reg.rimpl(X, y, alpha = alpha, lambda = lambda,
                                        init.coef = init.coef, maxit = ctrl$pense.maxit,
                                        control = ctrl)

    pr.c <- pense:::pen.s.reg(X, y, alpha = alpha, lambda = lambda,
                                  init.coef = init.coef, maxit = ctrl$pense.maxit,
                                  control = ctrl)

    expect_equal(pr.c$intercept, pr.r$intercept)
    expect_equal(pr.c$beta, pr.r$beta)
    expect_equal(pr.c$objF, pr.r$objF)


    ## EN penalty
    alpha <- 0.5

    pr.r <- pense:::pen.s.reg.rimpl(X, y, alpha = alpha, lambda = lambda,
                                        init.coef = init.coef, maxit = ctrl$pense.maxit,
                                        control = ctrl)

    pr.c <- pense:::pen.s.reg(X, y, alpha = alpha, lambda = lambda,
                                  init.coef = init.coef, maxit = ctrl$pense.maxit,
                                  control = ctrl)

    expect_equal(pr.c$intercept, pr.r$intercept)
    expect_equal(pr.c$beta, pr.r$beta)
    expect_equal(pr.c$objF, pr.r$objF)

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

    pr.r <- pense:::pen.s.reg.rimpl(X, y, alpha = alpha, lambda = lambda,
                                        init.coef = init.coef, maxit = ctrl$pense.maxit,
                                        control = ctrl)

    pr.c <- pense:::pen.s.reg(X, y, alpha = alpha, lambda = lambda,
                                  init.coef = init.coef, maxit = ctrl$pense.maxit,
                                  control = ctrl)

    expect_equal(pr.c$intercept, pr.r$intercept)
    expect_equal(pr.c$beta, pr.r$beta)
    expect_equal(pr.c$objF, pr.r$objF)


    ## EN penalty
    alpha <- 0.5

    pr.r <- pense:::pen.s.reg.rimpl(X, y, alpha = alpha, lambda = lambda,
                                        init.coef = init.coef, maxit = ctrl$pense.maxit,
                                        control = ctrl)

    pr.c <- pense:::pen.s.reg(X, y, alpha = alpha, lambda = lambda,
                                  init.coef = init.coef, maxit = ctrl$pense.maxit,
                                  control = ctrl)

    expect_equal(pr.c$intercept, pr.r$intercept)
    expect_equal(pr.c$beta, pr.r$beta)
    expect_equal(pr.c$objF, pr.r$objF)

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

    pr.r <- pense:::pen.s.reg.rimpl(X, y, alpha = alpha, lambda = lambda,
                                        init.coef = init.coef, maxit = ctrl$pense.maxit,
                                        control = ctrl)

    pr.c <- pense:::pen.s.reg(X, y, alpha = alpha, lambda = lambda,
                                  init.coef = init.coef, maxit = ctrl$pense.maxit,
                                  control = ctrl)

    expect_equal(pr.c$intercept, pr.r$intercept)
    expect_equal(pr.c$beta, pr.r$beta)
    expect_equal(pr.c$objF, pr.r$objF)


    ## EN penalty
    alpha <- 0.5

    pr.r <- pense:::pen.s.reg.rimpl(X, y, alpha = alpha, lambda = lambda,
                                        init.coef = init.coef, maxit = ctrl$pense.maxit,
                                        control = ctrl)

    pr.c <- pense:::pen.s.reg(X, y, alpha = alpha, lambda = lambda,
                                  init.coef = init.coef, maxit = ctrl$pense.maxit,
                                  control = ctrl)

    expect_equal(pr.c$intercept, pr.r$intercept)
    expect_equal(pr.c$beta, pr.r$beta)
    expect_equal(pr.c$objF, pr.r$objF)

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
        pense:::pen.s.reg.rimpl(X, y, alpha = alpha, lambda = lambda,
                                    init.coef = init.coef, maxit = ctrl$pense.maxit,
                                    control = ctrl),
        regexp = 'did not converge'
    )

    expect_warning(
        pense:::pen.s.reg(X, y, alpha = alpha, lambda = lambda,
                                      init.coef = init.coef, maxit = ctrl$pense.maxit,
                                      control = ctrl),
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

    pr.r <- pense:::pen.s.reg.rimpl(X, y, alpha = alpha, lambda = lambda,
                                        init.coef = init.coef, maxit = ctrl$pense.maxit,
                                        control = ctrl)

    pr.c <- pense:::pen.s.reg(X, y, alpha = alpha, lambda = lambda,
                                  init.coef = init.coef, maxit = ctrl$pense.maxit,
                                  control = ctrl)

    expect_equal(pr.c$intercept, pr.r$intercept)
    expect_equal(pr.c$beta, pr.r$beta)
    expect_equal(pr.c$objF, pr.r$objF)


    ## EN penalty
    alpha <- 0.5

    pr.r <- pense:::pen.s.reg.rimpl(X, y, alpha = alpha, lambda = lambda,
                                        init.coef = init.coef, maxit = ctrl$pense.maxit,
                                        control = ctrl)

    pr.c <- pense:::pen.s.reg(X, y, alpha = alpha, lambda = lambda,
                                  init.coef = init.coef, maxit = ctrl$pense.maxit,
                                  control = ctrl)

    expect_equal(pr.c$intercept, pr.r$intercept)
    expect_equal(pr.c$beta, pr.r$beta)
    expect_equal(pr.c$objF, pr.r$objF)

})
