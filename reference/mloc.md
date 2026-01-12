# Compute the M-estimate of Location

Compute the M-estimate of location using an auxiliary estimate of the
scale.

## Usage

``` r
mloc(x, scale, rho = "bisquare", eff = 0.9, cc, max_it = 200, eps = 1e-08)
```

## Arguments

- x:

  numeric values. Missing values are verbosely ignored.

- scale:

  scale of the `x` values. If omitted, uses the
  [mad()](https://rdrr.io/r/stats/mad.html).

- rho:

  the \\\rho\\ function to use. See [`rho_function()`](rho_function.md)
  for available functions.

- eff:

  desired efficiency under the Normal model.

- cc:

  value of the tuning constant for the chosen \\\rho\\ function. If
  specified, overrides the desired efficiency.

- max_it:

  maximum number of iterations.

- eps:

  numerical tolerance to check for convergence.

## Value

a single numeric value, the M-estimate of location.

## See also

Other functions to compute robust estimates of location and scale:
[`mlocscale()`](mlocscale.md), [`mscale()`](mscale.md),
[`tau_size()`](tau_size.md)
