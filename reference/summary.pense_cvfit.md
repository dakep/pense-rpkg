# Summarize Cross-Validated PENSE Fit

If `lambda = "se"` and `object` contains fitted estimates for every
penalization level in the sequence, extract the coefficients of the most
parsimonious model with prediction performance statistically
indistinguishable from the best model. This is determined to be the
model with prediction performance within `se_mult * cv_se` from the best
model.

## Usage

``` r
# S3 method for class 'pense_cvfit'
summary(object, alpha, lambda = "min", se_mult = 1, ...)

# S3 method for class 'pense_cvfit'
print(x, alpha, lambda = "min", se_mult = 1, ...)
```

## Arguments

- object, x:

  an (adaptive) PENSE fit with cross-validation information.

- alpha:

  Either a single number or missing. If given, only fits with the given
  `alpha` value are considered. If `lambda` is a numeric value and
  `object` was fit with multiple `alpha` values, the parameter `alpha`
  must not be missing.

- lambda:

  either a string specifying which penalty level to use (`"min"`,
  `"se"`, `"{x}-se`") or a single numeric value of the penalty
  parameter. See details.

- se_mult:

  If `lambda = "se"`, the multiple of standard errors to tolerate.

- ...:

  ignored.

## See also

[`prediction_performance()`](prediction_performance.md) for information
about the estimated prediction performance.

[`coef.pense_cvfit()`](coef.pense_cvfit.md) for extracting only the
estimated coefficients.

Other functions for plotting and printing:
[`plot.pense_cvfit()`](plot.pense_cvfit.md),
[`plot.pense_fit()`](plot.pense_fit.md),
[`prediction_performance()`](prediction_performance.md)
