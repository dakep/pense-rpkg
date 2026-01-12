# Extract Residuals

Extract residuals from a PENSE (or LS-EN) regularization path with
hyper-parameters chosen by cross-validation.

## Usage

``` r
# S3 method for class 'pense_cvfit'
residuals(object, alpha = NULL, lambda = "min", se_mult = 1, ...)
```

## Arguments

- object:

  PENSE with cross-validated hyper-parameters to extract coefficients
  from.

- alpha:

  Either a single number or `NULL` (default). If given, only fits with
  the given `alpha` value are considered. If `lambda` is a numeric value
  and `object` was fit with multiple *alpha* values and no value is
  provided, the first value in `object$alpha` is used with a warning.

- lambda:

  either a string specifying which penalty level to use (`"min"`,
  `"se"`, `"{m}-se`") or a single numeric value of the penalty
  parameter. See details.

- se_mult:

  If `lambda = "se"`, the multiple of standard errors to tolerate.

- ...:

  currently not used.

## Value

a numeric vector of residuals for the given penalization level.

## Hyper-parameters

If `lambda = "{m}-se"` and `object` contains fitted estimates for every
penalization level in the sequence, use the fit the most parsimonious
model with prediction performance statistically indistinguishable from
the best model. This is determined to be the model with prediction
performance within `m * cv_se` from the best model. If `lambda = "se"`,
the multiplier *m* is taken from `se_mult`.

By default all *alpha* hyper-parameters available in the fitted object
are considered. This can be overridden by supplying one or multiple
values in parameter `alpha`. For example, if `lambda = "1-se"` and
`alpha` contains two values, the "1-SE" rule is applied individually for
each `alpha` value, and the fit with the better prediction error is
considered.

In case `lambda` is a number and `object` was fit for several *alpha*
hyper-parameters, `alpha` must also be given, or the first value in
`object$alpha` is used with a warning.

## See also

Other functions for extracting components:
[`coef.pense_cvfit()`](coef.pense_cvfit.md),
[`coef.pense_fit()`](coef.pense_fit.md),
[`predict.pense_cvfit()`](predict.pense_cvfit.md),
[`predict.pense_fit()`](predict.pense_fit.md),
[`residuals.pense_fit()`](residuals.pense_fit.md)

## Examples

``` r
# Compute the LS-EN regularization path for Freeny's revenue data
# (see ?freeny)
data(freeny)
x <- as.matrix(freeny[ , 2:5])

regpath <- elnet(x, freeny$y, alpha = 0.75)

# Predict the response using a specific penalization level
predict(regpath, newdata = freeny[1:5, 2:5],
        lambda = regpath$lambda[[1]][[10]])
#>  1962.25   1962.5  1962.75     1963  1963.25 
#> 9.071638 9.075877 9.082341 9.091051 9.103643 

# Extract the residuals at a certain penalization level
residuals(regpath, lambda = regpath$lambda[[1]][[5]])
#>              Qtr1         Qtr2         Qtr3         Qtr4
#> 1962              -0.396169224 -0.398919454 -0.378220695
#> 1963 -0.384593892 -0.294338479 -0.278030430 -0.259140970
#> 1964 -0.265543684 -0.219857551 -0.206023748 -0.173545125
#> 1965 -0.192086388 -0.146253918 -0.133786817 -0.096671191
#> 1966 -0.090347926 -0.044601264 -0.027307275 -0.015531362
#> 1967  0.006739697  0.037235389  0.038869333  0.071407419
#> 1968  0.086275062  0.102036759  0.141643505  0.167310507
#> 1969  0.175437660  0.211500792  0.222942549  0.260029329
#> 1970  0.247870647  0.296181577  0.294470169  0.278217170
#> 1971  0.306686546  0.334684047  0.353729154  0.367702082

# Select penalization level via cross-validation
set.seed(123)
cv_results <- elnet_cv(x, freeny$y, alpha = 0.5,
                       cv_repl = 10, cv_k = 4)

# Predict the response using the "best" penalization level
predict(cv_results, newdata = freeny[1:5, 2:5])
#>  1962.25   1962.5  1962.75     1963  1963.25 
#> 8.795162 8.807070 8.824535 8.842175 8.882970 

# Extract the residuals at the "best" penalization level
residuals(cv_results)
#>              Qtr1         Qtr2         Qtr3         Qtr4
#> 1962              -0.002801588 -0.015699794 -0.009674689
#> 1963 -0.029165027  0.024539756  0.012794692  0.013400637
#> 1964 -0.014177641  0.006727561 -0.001786569  0.012847301
#> 1965 -0.026214283 -0.003738367 -0.010414370 -0.002497682
#> 1966 -0.014953173  0.014923561  0.008857587  0.002955549
#> 1967 -0.002438578  0.011336771 -0.004817002 -0.002203485
#> 1968 -0.014622332 -0.016038521  0.004783442  0.009509723
#> 1969  0.008935501  0.025135859  0.011008309  0.028171208
#> 1970 -0.021812388  0.015300404 -0.005936032 -0.018484409
#> 1971 -0.011338929  0.005850584  0.003783641  0.007952774
# Extract the residuals at a more parsimonious penalization level
residuals(cv_results, lambda = "1.5-se")
#>               Qtr1          Qtr2          Qtr3          Qtr4
#> 1962               -0.0133914763 -0.0253025812 -0.0180073467
#> 1963 -0.0375085950  0.0197904960  0.0063674275  0.0071084444
#> 1964 -0.0199445085  0.0028559476 -0.0059886764  0.0089958759
#> 1965 -0.0300906249 -0.0052788324 -0.0128420949 -0.0036364879
#> 1966 -0.0166186511  0.0137904162  0.0073750895  0.0016121587
#> 1967 -0.0021618673  0.0116385091 -0.0046227113  0.0002016235
#> 1968 -0.0117185912 -0.0127504673  0.0088610048  0.0131305182
#> 1969  0.0114094522  0.0285319822  0.0147126545  0.0325247699
#> 1970 -0.0158575338  0.0222103711  0.0004937030 -0.0127802705
#> 1971 -0.0035765461  0.0133307584  0.0117139290  0.0154227310
```
