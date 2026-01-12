# Changelog

## pense 2.5.0

- Use Robust Information Sharing Cross-Validation (RIS-CV) by default
  for [`pense_cv()`](../reference/pense_cv.md) and
  [`adapense_cv()`](../reference/pense_cv.md). Details about RIS-CV can
  be found in [Kepplinger & Wei
  (2025)](https://doi.org/10.1080/00401706.2025.2540970).
- Add support for the optimal rho function (Maronna et al. 2018, Section
  5.8.1), using the polynomial approximation from the robustbase
  package.
- Fix undefined behavior identified by the extended CRAN checks.

## pense 2.2.2

CRAN release: 2024-07-27

- Fix build warnings from GCC 14.1.0 on CRAN

## pense 2.2.1

- Fix build error on Windows using GCC version \< 9.

## pense 2.2.0

CRAN release: 2023-02-07

- Fix bug where argument `max_solutions` is not used correctly.
- Add new numerical coordinate-descent algorithms for LS-EN
  ([`en_cd_options()`](../reference/en_cd_options.md)).
- Revamp the regularization path and how starting points are used and
  shared.
- Share starting points across CV folds and replications.
- Address changes in the upcoming Rcpp release.

## pense 2.1.0

CRAN release: 2021-07-07

- Penalty loadings are now applied to both the L1 and L2 parts of the EN
  penalty. This will lead to different results for adaptive PENSE and
  other adaptive estimators when fitted with *alpha \< 1*!
- PENSE and regularized M-estimators now accept multiple `alpha` values
  and automatic hyper-parameter selection will also choose the best
  `alpha` value.
- New support for specifying the more general “1-SE” rules for the
  penalization level as string. All methods which support
  `lambda = "min"` to extract the best fit, also support the syntax
  `lambda = "{m}-se"` to extract the most parsimonious fit within *m*
  standard-errors of the best fit.
- Adaptively choose the actual breakdown point based on the number of
  observations. The chosen breakdown point is close to the
  user-specified breakdown point, but avoids numerical instabilities in
  the S-loss and excessive computation time caused by these
  instabilities.
- Simplify the DAL algorithm to fully rely on linear algebra routines
  from the BLAS/LAPACK library linked to R. To improve the speed of the
  DAL algorithm, optimized BLAS/LAPACK libraries are recommended.
- Fix memory issues from edge-cases and OpenMP problems with Intel
  compilers

## pense 2.0.3

CRAN release: 2021-04-14

- Fix a bug causing PENSE-Ridge, i.e., `pense(..., alpha = 0)`, to take
  a long time to compute.
- Fix a compilation error on RHEL due to an error in the autoconf
  script.
- Fix problems in
  [`prediction_performance()`](../reference/prediction_performance.md)
  related to the non-standard evaluation of objects.
- Also return standardized coefficients as `std_beta` and
  `std_intercept`. \# pense 2.0.2
- Fix mishandling of response variables with a robust scale of 0, e.g.,
  0-inflated responses or responses with more than 50% identical values.
  \# pense 2.0.1
- Add new functions for compute adaptive PENSE estimates
  ([`adapense()`](../reference/pense.md) and
  [`adapense_cv()`](../reference/pense_cv.md)).
- Functions for fitting the model ([`pense()`](../reference/pense.md),
  [`adapense()`](../reference/pense.md),
  [`regmest()`](../reference/regmest.md), etc.) are not estimating
  prediction performance via cross-validation anymore. This can now be
  done using the corresponding functions
  [`pense_cv()`](../reference/pense_cv.md),
  [`adapense_cv()`](../reference/pense_cv.md), and so on.
- New function
  [`prediction_performance()`](../reference/prediction_performance.md)
  to summarize the prediction performance of several fits.
- The [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`predict()`](https://rdrr.io/r/stats/predict.html) methods for
  cross-validated fits also implement the “one-standard-error rule”
  (with the “1” adjustable by the user).
- Decrease computation time for most problems.
- New ADMM algorithm for (weighted) elastic net problems with many
  observations and many predictors. The new algorithm can be selected
  with [`en_admm_options()`](../reference/en_admm_options.md).
- Argument `correct` in [`pense()`](../reference/pense.md), `pensem()`,
  [`coef()`](https://rdrr.io/r/stats/coef.html), etc., is not supported
  anymore and will be ignored with a warning. All estimates are now
  **uncorrected** (i.e., `correct=FALSE` in previous versions of the
  package).
- Make interface more consistent and deprecate the following methods:
  - `pensem()` is now called `pensem_cv()`.
  - `initest_options()` is replaced by
    [`enpy_options()`](../reference/enpy_options.md) using better naming
    of arguments.
  - `en_options_aug_lars()` and `en_options_dal()` are replaced by
    [`en_lars_options()`](../reference/en_lars_options.md) and
    [`en_dal_options()`](../reference/en_dal_options.md) for more
    consistent naming.
  - `pense_options()` and `mstep_options()` are superseded by
    [`mm_algorithm_options()`](../reference/mm_algorithm_options.md) and
    arguments specified in the calls to
    [`pense()`](../reference/pense.md) and companions.
  - `enpy()` is replaced by
    [`enpy_initial_estimates()`](../reference/enpy_initial_estimates.md)
    which has different default argument values.
- Deprecated functions can still be used (for now) with a warning.

## pense 1.2.9

CRAN release: 2020-02-09

- Fix LTO warnings reported in CRAN checks
- Update autoconf script to address deprecation warnings in r-devel.

## pense 1.2.5

CRAN release: 2019-06-08

- Fix compatibility of BLAS/LAPACK prototypes with RcppArmadillo
  0.9.500.

## pense 1.2.4

CRAN release: 2019-04-27

- Fix autoconf script.

## pense 1.2.1

CRAN release: 2019-01-17

- Prepare for changes to the upcoming *Rcpp* (make compatible with
  `STRICT_R_HEADERS`)
- Fix a bug in computing PSCs when using the augmented ridge algorithm
  for EN.

## pense 1.2.0

CRAN release: 2018-03-11

- Changed the internal scaling of the regularization parameter for
  `pense` and `pensem`. **Note**: The *lambda* values in this release
  are not the same as in previous releases!
- Fixed a bug when standardizing predictor variables with a MAD of 0
  (thanks [@hadjipantelis](https://github.com/hadjipantelis) for
  reporting).
- The maximum value for the regularization parameter lambda is now
  chosen exactly.
- Fixed a bug when computing “exact” principal sensitivity components.
  \# pense 1.0.8
- Fix error with robustbase-0.92-8 as reported by Martin Maechler.
- Fix undefined behavior in C++ code resulting in build error on Solaris
  (x86).
- Fix [`predict()`](https://rdrr.io/r/stats/predict.html) function for
  `pensem` objects if computed from a fitted `pense` object.
- Always use `delta` and `cc` specified in `pense_options()` for the
  initial estimator. Remove `delta` and `cc` arguments from
  `initest_options()` and instead add them to `enpy()`.
- Add further measure of the prediction performance (`resid_size`) to
  `obj$cv_lambda_grid`, where `obj` is of class `pense` or `pensem`. \#
  pense 1.0.6:
- Initial stable release of the package.
