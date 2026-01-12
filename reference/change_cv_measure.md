# Change the Cross-Validation Measure

For cross-validated fits using the RIS-CV strategy, the measure of
prediction accuracy can be adjusted post-hoc.

## Usage

``` r
change_cv_measure(
  x,
  measure = c("wrmspe", "wmape", "tau_size", "wrmspe_cv", "wmape_cv"),
  max_solutions = Inf
)
```

## Arguments

- x:

  fitted (adaptive) PENSE or M-estimator

- measure:

  the measure to use for prediction accuracy

- max_solutions:

  consider only this many of the best solutions. If missing, all
  solutions are considered.

## Value

a `pense.cvfit` object using the updated measure of prediction accuracy

## See also

Other functions to compute robust estimates with CV:
[`pense_cv()`](pense_cv.md), [`regmest_cv()`](regmest_cv.md)
