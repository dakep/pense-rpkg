# Get the constant for the desired efficiency of the M-estimate of location using the optimal \\\rho\\ function

Get the constant for the desired efficiency of the M-estimate of
location using the optimal \\\rho\\ function

## Usage

``` r
.mopt_efficiency_const(eff, eps = sqrt(.Machine$double.eps))
```

## Arguments

- eff:

  desired efficiency (between 0 and 1)

- eps:

  numerical tolerance for equality comparisons

## Value

tuning constant for desired efficiency
