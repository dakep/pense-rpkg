# ENPY Initial Estimates for EN S-Estimators

Compute initial estimates for the EN S-estimator using the EN-PY
procedure.

## Usage

``` r
enpy_initial_estimates(
  x,
  y,
  alpha,
  lambda,
  bdp = 0.25,
  cc,
  intercept = TRUE,
  penalty_loadings,
  enpy_opts = enpy_options(),
  mscale_opts = mscale_algorithm_options(),
  eps = 1e-06,
  sparse = FALSE,
  ncores = 1L
)
```

## Arguments

- x:

  `n` by `p` matrix of numeric predictors.

- y:

  vector of response values of length `n`.

- alpha:

  elastic net penalty mixing parameter with \\0 \le \alpha \le 1\\.
  `alpha = 1` is the LASSO penalty, and `alpha = 0` the Ridge penalty.
  Can be a vector of several values, but `alpha = 0` cannot be mixed
  with other values.

- lambda:

  a vector of positive values of penalization levels.

- bdp:

  desired breakdown point of the estimator, between 0.05 and 0.5. The
  actual breakdown point may be slightly larger/smaller to avoid
  instabilities of the S-loss.

- cc:

  cutoff value for the bisquare rho function. By default, chosen to
  yield a consistent estimate for the Normal distribution.

- intercept:

  include an intercept in the model.

- penalty_loadings:

  a vector of positive penalty loadings (a.k.a. weights) for different
  penalization of each coefficient. Only allowed for `alpha` \> 0.

- enpy_opts:

  options for the EN-PY algorithm, created with the
  [`enpy_options()`](enpy_options.md) function.

- mscale_opts:

  options for the M-scale estimation. See
  [`mscale_algorithm_options()`](mscale_algorithm_options.md) for
  details.

- eps:

  numerical tolerance.

- sparse:

  use sparse coefficient vectors.

- ncores:

  number of CPU cores to use in parallel. By default, only one CPU core
  is used. Not supported on all platforms, in which case a warning is
  given.

## Details

If these manually computed initial estimates are intended as starting
points for [`pense()`](pense.md), they are by default *shared* for all
penalization levels. To restrict the use of the initial estimates to the
penalty level they were computed for, use
`as_starting_point(..., specific = TRUE)`. See
[`as_starting_point()`](starting_point.md) for details.

## References

Cohen Freue, G.V.; Kepplinger, D.; Salibián-Barrera, M.; Smucler, E.
Robust elastic net estimators for variable selection and identification
of proteomic biomarkers. *Ann. Appl. Stat.* **13** (2019), no. 4,
2065–2090 [doi:10.1214/19-AOAS1269](https://doi.org/10.1214/19-AOAS1269)

## See also

Other functions for initial estimates:
[`enpy_options()`](enpy_options.md), [`prinsens()`](prinsens.md),
[`starting_point()`](starting_point.md)
