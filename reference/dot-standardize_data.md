# Standardize data

Standardize data

## Usage

``` r
.standardize_data(
  x,
  y,
  intercept,
  standardize,
  robust,
  sparse,
  mscale_opts,
  cc,
  target_scale_x = NULL,
  ...
)
```

## Arguments

- x:

  predictor matrix. Can also be a list with components `x` and `y`, in
  which case `y` is ignored.

- y:

  response vector.

- intercept:

  is an intercept included (i.e., should `y` be centered?)

- standardize:

  standardize or not.

- robust:

  use robust standardization.

- cc:

  cutoff value for the rho functions used in scale and location
  estimates.

- ...:

  passed on to [`mlocscale()`](mlocscale.md).

## Value

a list with the following entries:
