# Compute the M-Scale of Centered Values

Compute the M-scale without centering the values.

## Usage

``` r
mscale(x, bdp = 0.25, cc, opts = mscale_algorithm_options())
```

## Arguments

- x:

  numeric values. Missing values are verbosely ignored.

- bdp:

  desired breakdown point (between 0 and 0.5).

- cc:

  tuning parameters for the chosen rho function. By default, chosen to
  yield a consistent estimate for the Normal distribution.

- opts:

  a list of options for the M-scale estimation algorithm, see
  [`mscale_algorithm_options()`](mscale_algorithm_options.md) for
  details.

## Value

the M-estimate of scale.

## See also

Other functions to compute robust estimates of location and scale:
[`mloc()`](mloc.md), [`mlocscale()`](mlocscale.md),
[`tau_size()`](tau_size.md)
