# Create Starting Points for the PENSE Algorithm

Create a starting point for starting the PENSE algorithm in
[`pense()`](pense.md). Multiple starting points can be created by
combining starting points via
`c(starting_point_1, starting_point_2, ...)`.

## Usage

``` r
starting_point(beta, intercept, lambda, alpha)

as_starting_point(object, specific = FALSE, ...)

# S3 method for class 'enpy_starting_points'
as_starting_point(object, specific = FALSE, ...)

# S3 method for class 'pense_fit'
as_starting_point(object, specific = FALSE, alpha, lambda, ...)

# S3 method for class 'pense_cvfit'
as_starting_point(
  object,
  specific = FALSE,
  alpha,
  lambda = c("min", "se"),
  se_mult = 1,
  ...
)
```

## Arguments

- beta:

  beta coefficients at the starting point. Can be a numeric vector, a
  sparse vector of class
  [dsparseVector](https://rdrr.io/pkg/Matrix/man/sparseVector-class.html),
  or a sparse matrix of class
  [dgCMatrix](https://rdrr.io/pkg/Matrix/man/CsparseMatrix-class.html)
  with a single column.

- intercept:

  intercept coefficient at the starting point.

- lambda:

  optionally either a string specifying which penalty level to use
  (`"min"` or `"se"`) or a numeric vector of the penalty levels to
  extract from `object`. Penalization levels not present in `object` are
  ignored with a warning. If `NULL`, all estimates in `object` are
  extracted. If a numeric vector, `alpha` must be given and a single
  number.

- alpha:

  optional value for the `alpha` hyper-parameter. If given, only
  estimates with matching `alpha` values are extracted. Values not
  present in `object` are ignored with a warning.

- object:

  an object with estimates to use as starting points.

- specific:

  whether the estimates should be used as starting points only at the
  penalization level they are computed for. Defaults to using the
  estimates as starting points for all penalization levels.

- ...:

  further arguments passed to or from other methods.

- se_mult:

  If `lambda = "se"`, the multiple of standard errors to tolerate.

## Value

an object of type `starting_points` to be used as starting point for
[`pense()`](pense.md).

## Details

A starting points can either be *shared*, i.e., used for every
penalization level PENSE estimates are computed for, or *specific* to
one penalization level. To create a specific starting point, provide the
penalization parameters `lambda` and `alpha`. If `lambda` or `alpha` are
missing, a shared starting point is created. Shared and specific
starting points can all be combined into a single list of starting
points, with [`pense()`](pense.md) handling them correctly. Note that
specific starting points will lead to the `lambda` value being added to
the grid of penalization levels. See [`pense()`](pense.md) for details.

Starting points computed via
[`enpy_initial_estimates()`](enpy_initial_estimates.md) are by default
*shared* starting points but can be transformed to *specific* starting
points via `as_starting_point(..., specific = TRUE)`.

When creating starting points from cross-validated fits, it is possible
to extract only the estimate with best CV performance
(`lambda = "min"`), or the estimate with CV performance statistically
indistinguishable from the best performance (`lambda = "se"`). This is
determined to be the estimate with prediction performance within
`se_mult * cv_se` from the best model.

## See also

Other functions for initial estimates:
[`enpy_initial_estimates()`](enpy_initial_estimates.md),
[`enpy_options()`](enpy_options.md), [`prinsens()`](prinsens.md)
