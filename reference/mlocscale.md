# Compute the M-estimate of Location and Scale

Simultaneous estimation of the location and scale by means of
M-estimates.

## Usage

``` r
mlocscale(
  x,
  bdp = 0.25,
  eff = 0.9,
  scale_cc,
  location_rho,
  location_cc,
  opts = mscale_algorithm_options()
)
```

## Arguments

- x:

  numeric values. Missing values are verbosely ignored.

- bdp:

  desired breakdown point (between 0 and 0.5).

- eff:

  desired efficiency of the location estimate (between 0.1 and 0.99).

- scale_cc:

  tuning constant for the \\\rho\\ function for computing the scale
  estimate. By default, chosen to yield a consistent estimate for
  normally distributed values.

- location_rho:

  \\\rho\\ function for computing the location estimate. If missing, use
  the same function as for the scale estimate (`opts$rho`). See
  [`rho_function()`](rho_function.md) for a list of available \\\rho\\
  functions.

- location_cc:

  tuning constant for the location \\\rho\\ function. By default chosen
  to yield the desired efficiency. If this is provided, the desired
  efficiency is ignored.

- opts:

  a list of options for the M-scale estimating equations, See
  [`mscale_algorithm_options()`](mscale_algorithm_options.md) for
  details.

## Value

a vector with 2 elements, the M-estimate of location and the M-scale
estimate.

## See also

Other functions to compute robust estimates of location and scale:
[`mloc()`](mloc.md), [`mscale()`](mscale.md),
[`tau_size()`](tau_size.md)
