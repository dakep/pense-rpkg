# Cross-validation for (Adaptive) PENSE Estimates

Perform (repeated) K-fold cross-validation for [`pense()`](pense.md).

`adapense_cv()` is a convenience wrapper to compute adaptive PENSE
estimates.

## Usage

``` r
pense_cv(
  x,
  y,
  standardize = TRUE,
  lambda,
  cv_k,
  cv_repl = 1,
  cv_type = c("ris", "naive"),
  cv_metric = c("tau_size", "mape", "rmspe", "auroc"),
  ris_min_similarity = 0.5,
  fit_all = TRUE,
  fold_starts = c("full", "enpy", "both"),
  cv_algorithm_opts,
  cl = NULL,
  ...
)

adapense_cv(x, y, alpha, alpha_preliminary = 0, exponent = 1, ...)
```

## Arguments

- x:

  `n` by `p` matrix of numeric predictors.

- y:

  vector of response values of length `n`. For binary classification,
  `y` should be a factor with 2 levels.

- standardize:

  whether to standardize the `x` variables prior to fitting the PENSE
  estimates. Can also be set to `"cv_only"`, in which case the input
  data is not standardized, but the training data in the CV folds is
  scaled to match the scaling of the input data. Coefficients are always
  returned on the original scale. This can fail for variables with a
  large proportion of a single value (e.g., zero-inflated data). In this
  case, either compute with `standardize = FALSE` or standardize the
  data manually.

- lambda:

  optional user-supplied sequence of penalization levels. If given and
  not `NULL`, `nlambda` and `lambda_min_ratio` are ignored.

- cv_k:

  number of folds per cross-validation.

- cv_repl:

  number of cross-validation replications.

- cv_type:

  what kind of cross-validation should be performed: robust information
  sharing (`ris`) or standard (`naive`) CV.

- cv_metric:

  only for `cv_type='naive'`. Either a string specifying the performance
  metric to use, or a function to evaluate prediction errors in a single
  CV replication. If a function, the number of arguments define the data
  the function receives. If the function takes a single argument, it is
  called with a single numeric vector of prediction errors. If the
  function takes two or more arguments, it is called with the predicted
  values as first argument and the true values as second argument. The
  function must always return a single numeric value quantifying the
  prediction performance. The order of the given values corresponds to
  the order in the input data.

- ris_min_similarity:

  minimum average similarity of the CV solutions to be considered
  (between 0 and 1). If no CV solution satisfies this lower bound, the
  best CV solution will be used regardless of similarity.

- fit_all:

  only for `cv_type='naive'`. If `TRUE`, fit the model for all
  penalization levels. Can also be any combination of `"min"` and
  `"{x}-se"`, in which case only models at the penalization level with
  smallest average CV accuracy, or within `{x}` standard errors,
  respectively. Setting `fit_all` to `FALSE` is equivalent to `"min"`.
  Applies to all `alpha` value.

- fold_starts:

  how to determine starting values in the cross-validation folds. If
  `"full"` (default), use the best solution from the fit to the full
  data as starting value. This implies `fit_all=TRUE`. If `"enpy"`
  compute separate ENPY initial estimates in each fold. The option
  `"both"` uses both. These starts are in addition to the starts
  provided in `other_starts`.

- cv_algorithm_opts:

  Override algorithm options for the CV iterations. This is usually not
  necessary, unless the user wants to change the number of solutions
  retained for the CV training data.

- cl:

  a [parallel](https://rdrr.io/r/parallel/makeCluster.html) cluster. Can
  only be used in combination with `ncores = 1`.

- ...:

  Arguments passed on to [`pense`](pense.md)

  `nlambda`

  :   number of penalization levels.

  `lambda_min_ratio`

  :   Smallest value of the penalization level as a fraction of the
      largest level (i.e., the smallest value for which all coefficients
      are zero). The default depends on the sample size relative to the
      number of variables and `alpha`. If more observations than
      variables are available, the default is `1e-3 * alpha`, otherwise
      `1e-2 * alpha`.

  `nlambda_enpy`

  :   number of penalization levels where the EN-PY initial estimate is
      computed.

  `penalty_loadings`

  :   a vector of positive penalty loadings (a.k.a. weights) for
      different penalization of each coefficient. Only allowed for
      `alpha` \> 0.

  `enpy_lambda`

  :   optional user-supplied sequence of penalization levels at which
      EN-PY initial estimates are computed. If given and not `NULL`,
      `nlambda_enpy` is ignored.

  `other_starts`

  :   a list of other staring points, created by
      [`starting_point()`](starting_point.md). If the output of
      [`enpy_initial_estimates()`](enpy_initial_estimates.md) is given,
      the starting points will be *shared* among all penalization
      levels. Note that if a the starting point is *specific* to a
      penalization level, this penalization level is added to the grid
      of penalization levels (either the manually specified grid in
      `lambda` or the automatically generated grid of size `nlambda`).
      If `standardize = TRUE`, the starting points are also scaled.

  `intercept`

  :   include an intercept in the model.

  `bdp`

  :   desired breakdown point of the estimator, between 0.05 and 0.5.
      The actual breakdown point may be slightly larger/smaller to avoid
      instabilities of the S-loss.

  `cc`

  :   tuning constant for the S-estimator. Default is chosen based on
      the breakdown point `bdp`. This affects the estimated coefficients
      only if `standardize=TRUE`. Otherwise only the estimated scale of
      the residuals would be affected.

  `eps`

  :   numerical tolerance.

  `explore_solutions`

  :   number of solutions to keep after the exploration step. The best
      `explore_solutions` are then iterated to full numerical tolerance
      `eps`. If 0, all non-duplicated solutions are kept.

  `explore_tol,explore_it`

  :   numerical tolerance and maximum number of iterations for exploring
      possible solutions. The tolerance should be (much) looser than
      `eps` to be useful, and the number of iterations should also be
      much smaller than the maximum number of iterations given via
      `algorithm_opts`. `explore_tol` is also used to determine if two
      solutions are equal in the exploration stage.

  `max_solutions`

  :   retain only up to `max_solutions` unique solutions per
      penalization level.

  `comparison_tol`

  :   numeric tolerance to determine if two solutions are equal. The
      comparison is first done on the absolute difference in the value
      of the objective function at the solution. If this is less than
      `comparison_tol`, two solutions are deemed equal if the squared
      difference of the intercepts is less than `comparison_tol` and the
      squared \\L_2\\ norm of the difference vector is less than
      `comparison_tol`.

  `add_zero_based`

  :   also consider the 0-based regularization path. See details for a
      description.

  `enpy_specific`

  :   use the EN-PY initial estimates only at the penalization level
      they are computed for. See details for a description.

  `carry_forward`

  :   carry the best solutions forward to the next penalty level.

  `sparse`

  :   use sparse coefficient vectors.

  `ncores`

  :   number of CPU cores to use in parallel. By default, only one CPU
      core is used. Not supported on all platforms, in which case a
      warning is given.

  `algorithm_opts`

  :   options for the MM algorithm to compute the estimates. See
      [`mm_algorithm_options()`](mm_algorithm_options.md) for details.

  `mscale_opts`

  :   options for the M-scale estimation. See
      [`mscale_algorithm_options()`](mscale_algorithm_options.md) for
      details.

  `enpy_opts`

  :   options for the ENPY initial estimates, created with the
      [`enpy_options()`](enpy_options.md) function. See
      [`enpy_initial_estimates()`](enpy_initial_estimates.md) for
      details.

  `cv_k,cv_objective`

  :   deprecated and ignored. See `pense_cv()` for estimating prediction
      performance via cross-validation.

- alpha:

  elastic net penalty mixing parameter with \\0 \le \alpha \le 1\\.
  `alpha = 1` is the LASSO penalty, and `alpha = 0` the Ridge penalty.
  Can be a vector of several values, but `alpha = 0` cannot be mixed
  with other values.

- alpha_preliminary:

  `alpha` parameter for the preliminary estimate.

- exponent:

  the exponent for computing the penalty loadings based on the
  preliminary estimate.

## Value

a list-like object with the same components as returned by
[`pense()`](pense.md), plus the following:

- `cvres`:

  data frame of average cross-validated performance.

a list-like object as returned by `pense_cv()` plus the following

- `preliminary`:

  the CV results for the preliminary estimate.

- `exponent`:

  exponent used to compute the penalty loadings.

- `penalty_loadings`:

  penalty loadings used for the adaptive PENSE estimate.

## Details

The built-in CV metrics are

- `"tau_size"`:

  \\\tau\\-size of the prediction error, computed by
  [`tau_size()`](tau_size.md) (default).

- `"mape"`:

  Median absolute prediction error.

- `"rmspe"`:

  Root mean squared prediction error.

- `"auroc"`:

  Area under the receiver operator characteristic curve (actually 1 -
  AUROC). Only sensible for binary responses.

`adapense_cv()` is a convenience wrapper which performs 3 steps:

1.  compute preliminary estimates via
    `pense_cv(..., alpha = alpha_preliminary)`,

2.  computes the penalty loadings from the estimate `beta` with best
    prediction performance by
    `adapense_loadings = 1 / abs(beta)^exponent`, and

3.  compute the adaptive PENSE estimates via
    `pense_cv(..., penalty_loadings = adapense_loadings)`.

## See also

[`pense()`](pense.md) for computing regularized S-estimates without
cross-validation.

[`coef.pense_cvfit()`](coef.pense_cvfit.md) for extracting coefficient
estimates.

[`plot.pense_cvfit()`](plot.pense_cvfit.md) for plotting the CV
performance or the regularization path.

Other functions to compute robust estimates with CV:
[`change_cv_measure()`](change_cv_measure.md),
[`regmest_cv()`](regmest_cv.md)

Other functions to compute robust estimates with CV:
[`change_cv_measure()`](change_cv_measure.md),
[`regmest_cv()`](regmest_cv.md)

## Examples

``` r
# Compute the adaptive PENSE regularization path for Freeny's
# revenue data (see ?freeny)
data(freeny)
x <- as.matrix(freeny[ , 2:5])

## Either use the convenience function directly ...
set.seed(123)
ada_convenience <- adapense_cv(x, freeny$y, alpha = 0.5,
                               cv_repl = 2, cv_k = 4)

## ... or compute the steps manually:
# Step 1: Compute preliminary estimates with CV
set.seed(123)
preliminary_estimate <- pense_cv(x, freeny$y, alpha = 0,
                                 cv_repl = 2, cv_k = 4)
plot(preliminary_estimate, se_mult = 1)

# Step 2: Use the coefficients with best prediction performance
# to define the penalty loadings:
prelim_coefs <- coef(preliminary_estimate, lambda = 'min')
pen_loadings <- 1 / abs(prelim_coefs[-1])

# Step 3: Compute the adaptive PENSE estimates and estimate
# their prediction performance.
set.seed(123)
ada_manual <- pense_cv(x, freeny$y, alpha = 0.5,
                       cv_repl = 2, cv_k = 4,
                       penalty_loadings = pen_loadings)

# Visualize the prediction performance and coefficient path of
# the adaptive PENSE estimates (manual vs. automatic)
def.par <- par(no.readonly = TRUE)
layout(matrix(1:4, ncol = 2, byrow = TRUE))
plot(ada_convenience$preliminary)
plot(preliminary_estimate)
plot(ada_convenience)
plot(ada_manual)

par(def.par)
```
