# Principal Sensitivity Components

Compute Principal Sensitivity Components for Elastic Net Regression

## Usage

``` r
prinsens(
  x,
  y,
  alpha,
  lambda,
  intercept = TRUE,
  penalty_loadings,
  en_algorithm_opts,
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

  optional user-supplied sequence of penalization levels. If given and
  not `NULL`, `nlambda` and `lambda_min_ratio` are ignored.

- intercept:

  include an intercept in the model.

- penalty_loadings:

  a vector of positive penalty loadings (a.k.a. weights) for different
  penalization of each coefficient. Only allowed for `alpha` \> 0.

- en_algorithm_opts:

  options for the LS-EN algorithm. See
  [en_algorithm_options](en_algorithm_options.md) for details.

- eps:

  numerical tolerance.

- sparse:

  use sparse coefficient vectors.

- ncores:

  number of CPU cores to use in parallel. By default, only one CPU core
  is used. Not supported on all platforms, in which case a warning is
  given.

## Value

a list of principal sensitivity components, one per element in `lambda`.
Each PSC is itself a list with items `lambda`, `alpha`, and `pscs`.

## References

Cohen Freue, G.V.; Kepplinger, D.; Salibián-Barrera, M.; Smucler, E.
Robust elastic net estimators for variable selection and identification
of proteomic biomarkers. *Ann. Appl. Stat.* **13** (2019), no. 4,
2065–2090 [doi:10.1214/19-AOAS1269](https://doi.org/10.1214/19-AOAS1269)

Pena, D., and Yohai, V.J. A Fast Procedure for Outlier Diagnostics in
Large Regression Problems. *J. Amer. Statist. Assoc.* **94** (1999). no.
446, 434–445. [doi:10.2307/2670164](https://doi.org/10.2307/2670164)

## See also

Other functions for initial estimates:
[`enpy_initial_estimates()`](enpy_initial_estimates.md),
[`enpy_options()`](enpy_options.md),
[`starting_point()`](starting_point.md)
