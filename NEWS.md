# pense 2.0.2
  * Fix mishandling of response variables with a robust scale of 0, e.g., 0-inflated responses or responses with more than 50% identical values.
# pense 2.0.1
  * Add new functions for compute adaptive PENSE estimates (`adapense()` and `adapense_cv()`).
  * Functions for fitting the model (`pense()`, `adapense()`, `regmest()`, etc.) are not estimating prediction performance via cross-validation anymore.
    This can now be done using the corresponding functions `pense_cv()`, `adapense_cv()`, and so on.
  * New function `prediction_performance()` to summarize the prediction performance of several fits.
  * The `plot()`, `coef()`, `summary()`, and `predict()` methods for cross-validated fits also implement the "one-standard-error rule" (with the "1" adjustable by the user).
  * Decrease computation time for most problems.
  * New ADMM algorithm for (weighted) elastic net problems with many observations and many predictors.
    The new algorithm can be selected with `en_admm_options()`.
  * Argument `correct` in `pense()`, `pensem()`, `coef()`, etc., is not supported anymore and will be ignored with a warning.
    All estimates are now **uncorrected** (i.e., `correct=FALSE` in previous versions of the package).
  * Make interface more consistent and deprecate the following methods:
    - `pensem()` is now called `pensem_cv()`.
    - `initest_options()` is replaced by `enpy_options()` using better naming of arguments.
    - `en_options_aug_lars()` and `en_options_dal()` are replaced by `en_lars_options()` and `en_dal_options()` for more consistent naming.
    - `pense_options()` and `mstep_options()` are superseded by `mm_algorithm_options()` and arguments specified in the calls to `pense()` and companions.
    - `enpy()` is replaced by `enpy_initial_estimates()` which has different default argument values.
  * Deprecated functions can still be used (for now) with a warning.

# pense 1.2.9
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
