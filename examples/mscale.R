## Estimate the M-scale of a vector of values
set.seed(1234)
x <- rnorm(100)
mscale(x)
mscale(x, b = 0.25) # For a breakdown point of 25%
mscale(x, rho = "huber") # Using Huber's rho function
