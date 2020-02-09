# pense 1.2.8
  * Fix LTO warnings reported in CRAN checks
  * Update autoconf script to address deprecation warnings in r-devel.

# pense 1.2.5
  * Fix compatibility of BLAS/LAPACK prototypes with RcppArmadillo 0.9.500.

# pense 1.2.4
  * Fix autoconf script.

# pense 1.2.1
  * Prepare for changes to the upcoming _Rcpp_ (make compatible with `STRICT_R_HEADERS`)
  * Fix a bug in computing PSCs when using the augmented ridge algorithm for EN.

# pense 1.2.0
  * Changed the internal scaling of the regularization parameter for `pense` and `pensem`.
    **Note**: The _lambda_ values in this release are not the same as in previous releases!
  * Fixed a bug when standardizing predictor variables with a MAD of 0 (thanks @hadjipantelis for reporting).
  * The maximum value for the regularization parameter lambda is now chosen exactly.
  * Fixed a bug when computing "exact" principal sensitivity components.
# pense 1.0.8
  * Fix error with robustbase-0.92-8 as reported by Martin Maechler.
  * Fix undefined behavior in C++ code resulting in build error on Solaris (x86).
  * Fix `predict()` function for `pensem` objects if computed from a fitted `pense` object.
  * Always use `delta` and `cc` specified in `pense_options()` for the initial estimator. Remove `delta` and `cc` arguments from `initest_options()` and instead add them to `enpy()`.
  * Add further measure of the prediction performance (`resid_size`) to `obj$cv_lambda_grid`, where `obj` is of class `pense` or `pensem`.
# pense 1.0.6:
  * Initial stable release of the package.
