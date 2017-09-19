##
## A very simple example on artificial data
##

# Generate some dummy data
set.seed(1234)
n <- 50
p <- 40
x <- 1 + matrix(rnorm(n * p), ncol = p)
y <- x %*% c(2:5, numeric(p - 4)) + rnorm(n)

x_test <- 1 + matrix(rnorm(10 * n * p), ncol = p)
y_test <- x_test %*% c(2:5, numeric(p - 4)) + rnorm(n)

# Compute the S-estimator with an EN penalty for 30 lambda values
# (Note: In real applications, warm_reset should be at least 5)
set.seed(1234)
est <- pense(
    x, y,
    alpha = 0.8,
    nlambda = 30,
    warm_reset = 1L,
    cv_k = 5
)

# We can plot the CV prediction error curve
with(est,
     plot(
         cvavg ~ lambda,
         cv_lambda_grid,
         log = "x",
         col = 1L + (lambda == lambda_opt),
         pch = 20L,
         xlab = expression(lambda),
         ylab = "Robust RMSPE",
         main = "CV Prediction Error"
     )
)

# What is the RMSPE on test data
(rmspe <- sqrt(mean((y_test - predict(est, newdata = x_test))^2)))

##
## What happens if we replace 5 observations in the dummy data
## with outliers?
##
y[1:3] <- rnorm(3, -1000)

# Compute the S-estimator again
# (Note: In real applications, warm_reset should be at least 5)
set.seed(12345)
est_out <- pense(
    x, y,
    alpha = 0.8,
    nlambda = 30,
    warm_reset = 1L,
    cv_k = 5
)

# How does the RMSPE compare?
rmspe_out <- sqrt(mean((y_test - predict(est_out, newdata = x_test))^2))
c(rmspe = rmspe, rmspe_out = rmspe_out)

