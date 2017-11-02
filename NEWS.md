# pense 1.0.8
    * Fix bug reported by Martin Maechler when determining the largest lambda in the grid.
    * Fix undefined behaviour in C++ code resulting in build error on Solaris.
    * Fix `predict()` function for `pensem` objects if trained from a fitted `pense` object.
    * Always use `delta` and `cc` specified in `pense_options()` for the initial estimator. Remove `delta` and `cc` arguments from `initest_options()` and instead add them to `enpy()`.
    * Add further measure of the prediction performance (`resid_size`) to `obj$cv_lambda_grid`, where `obj` is of class `pense` or `pensem`.
# pense 1.0.6:
    * Initial stable release of the package.
