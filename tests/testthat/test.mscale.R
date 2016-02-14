test_that("M-scale", {
    expect_identical(penseinit::mscale(numeric(20)), NA_real_)
    set.seed(12345)
    expect_gt(penseinit::mscale(rnorm(2)), 0)
    expect_equal(penseinit::mscale(1 + numeric(10)),
                 penseinit::mscale(-1 + numeric(10))
    )
})
