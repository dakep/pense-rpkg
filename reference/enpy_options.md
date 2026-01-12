# Options for the ENPY Algorithm

Additional control options for the elastic net Peña-Yohai procedure.

## Usage

``` r
enpy_options(
  max_it = 10,
  keep_psc_proportion = 0.5,
  en_algorithm_opts,
  keep_residuals_measure = c("threshold", "proportion"),
  keep_residuals_proportion = 0.5,
  keep_residuals_threshold = 2,
  retain_best_factor = 2,
  retain_max = 500
)
```

## Arguments

- max_it:

  maximum number of EN-PY iterations.

- keep_psc_proportion:

  how many observations should to keep based on the Principal
  Sensitivity Components.

- en_algorithm_opts:

  options for the LS-EN algorithm. See
  [en_algorithm_options](en_algorithm_options.md) for details.

- keep_residuals_measure:

  how to determine what observations to keep, based on their residuals.
  If `proportion`, a fixed number of observations is kept. If
  `threshold`, only observations with residuals below the threshold are
  kept.

- keep_residuals_proportion:

  proportion of observations to kept based on their residuals.

- keep_residuals_threshold:

  only observations with (standardized) residuals less than this
  threshold are kept.

- retain_best_factor:

  only keep candidates that are within this factor of the best
  candidate. If `<= 1`, only keep candidates from the last iteration.

- retain_max:

  maximum number of candidates, i.e., only the best `retain_max`
  candidates are retained.

## Value

options for the ENPY algorithm.

## Details

The EN-PY procedure for computing initial estimates iteratively cleans
the data of observations with possibly outlying residual or high
leverage. Least-squares elastic net (LS-EN) estimates are computed on
the possibly clean subsets. At each iteration, the Principal Sensitivity
Components are computed to remove observations with potentially high
leverage. Among all the LS-EN estimates, the estimate with smallest
M-scale of the residuals is selected. Observations with largest residual
for the selected estimate are removed and the next iteration is started.

## See also

Other functions for initial estimates:
[`enpy_initial_estimates()`](enpy_initial_estimates.md),
[`prinsens()`](prinsens.md), [`starting_point()`](starting_point.md)
