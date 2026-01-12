# Prediction Performance of Adaptive PENSE Fits

Extract the prediction performance of one or more (adaptive) PENSE fits.

## Usage

``` r
prediction_performance(..., alpha = NULL, lambda = "min", se_mult = 1)

# S3 method for class 'pense_pred_perf'
print(x, ...)
```

## Arguments

- ...:

  one or more (adaptive) PENSE fits with cross-validation information.

- alpha:

  Either a numeric vector or `NULL` (default). If given, only fits with
  the given `alpha` value are considered. If `lambda` is a numeric value
  and `object` was fit with multiple `alpha` values, the parameter
  `alpha` must not be missing.

- lambda:

  either a string specifying which penalty level to use (`"min"`,
  `"se"`, `"{x}-se`") or a single numeric value of the penalty
  parameter. See details.

- se_mult:

  If `lambda = "se"`, the multiple of standard errors to tolerate.

- x:

  an object with information on prediction performance created with
  `prediction_performance()`.

## Value

a data frame with details about the prediction performance of the given
PENSE fits. The data frame has a custom print method summarizing the
prediction performances.

## Details

If `lambda = "se"` and the cross-validation was performed with multiple
replications, use the penalty level whit prediction performance within
`se_mult` of the best prediction performance.

## See also

[`summary.pense_cvfit()`](summary.pense_cvfit.md) for a summary of the
fitted model.

Other functions for plotting and printing:
[`plot.pense_cvfit()`](plot.pense_cvfit.md),
[`plot.pense_fit()`](plot.pense_fit.md),
[`summary.pense_cvfit()`](summary.pense_cvfit.md)
