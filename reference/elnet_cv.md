# Cross-validation for Least-Squares (Adaptive) Elastic Net Estimates

Perform (repeated) K-fold cross-validation for [`elnet()`](elnet.md).

## Usage

``` r
elnet_cv(
  x,
  y,
  lambda,
  cv_k,
  cv_repl = 1,
  cv_type = "naive",
  cv_metric = c("rmspe", "tau_size", "mape", "auroc"),
  fit_all = TRUE,
  cl = NULL,
  ...
)
```

## Arguments

- x:

  `n` by `p` matrix of numeric predictors.

- y:

  vector of response values of length `n`. For binary classification,
  `y` should be a factor with 2 levels.

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

- fit_all:

  only for `cv_type='naive'`. If `TRUE`, fit the model for all
  penalization levels. Can also be any combination of `"min"` and
  `"{x}-se"`, in which case only models at the penalization level with
  smallest average CV accuracy, or within `{x}` standard errors,
  respectively. Setting `fit_all` to `FALSE` is equivalent to `"min"`.
  Applies to all `alpha` value.

- cl:

  a [parallel](https://rdrr.io/r/parallel/makeCluster.html) cluster. Can
  only be used in combination with `ncores = 1`.

- ...:

  Arguments passed on to [`elnet`](elnet.md)

  `alpha`

  :   elastic net penalty mixing parameter with \\0 \le \alpha \le 1\\.
      `alpha = 1` is the LASSO penalty, and `alpha = 0` the Ridge
      penalty. Can be a vector of several values, but `alpha = 0` cannot
      be mixed with other values.

  `nlambda`

  :   number of penalization levels.

  `lambda_min_ratio`

  :   Smallest value of the penalization level as a fraction of the
      largest level (i.e., the smallest value for which all coefficients
      are zero). The default depends on the sample size relative to the
      number of variables and `alpha`. If more observations than
      variables are available, the default is `1e-3 * alpha`, otherwise
      `1e-2 * alpha`.

  `penalty_loadings`

  :   a vector of positive penalty loadings (a.k.a. weights) for
      different penalization of each coefficient.

  `standardize`

  :   standardize variables to have unit variance. Coefficients are
      always returned in original scale.

  `weights`

  :   a vector of positive observation weights.

  `intercept`

  :   include an intercept in the model.

  `sparse`

  :   use sparse coefficient vectors.

  `en_algorithm_opts`

  :   options for the EN algorithm. See
      [en_algorithm_options](en_algorithm_options.md) for details.

  `eps`

  :   numerical tolerance.

## Value

a list-like object with the same components as returned by
[`elnet()`](elnet.md), plus the following:

- `cvres`:

  data frame of average cross-validated performance.

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

## See also

[`elnet()`](elnet.md) for computing the LS-EN regularization path
without cross-validation.

[`pense_cv()`](pense_cv.md) for cross-validation of S-estimates of
regression with elastic net penalty.

[`coef.pense_cvfit()`](coef.pense_cvfit.md) for extracting coefficient
estimates.

[`plot.pense_cvfit()`](plot.pense_cvfit.md) for plotting the CV
performance or the regularization path.

Other functions for computing non-robust estimates:
[`elnet()`](elnet.md)

## Examples

``` r
# Compute the LS-EN regularization path for Freeny's revenue data
# (see ?freeny)
data(freeny)
x <- as.matrix(freeny[ , 2:5])

regpath <- elnet(x, freeny$y, alpha = c(0.5, 0.75))
plot(regpath)

plot(regpath, alpha = 0.75)


# Extract the coefficients at a certain penalization level
coef(regpath, lambda = regpath$lambda[[1]][[5]],
     alpha = 0.75)
#>           (Intercept) lag.quarterly.revenue           price.index 
#>              9.306304              0.000000              0.000000 
#>          income.level      market.potential 
#>              0.000000              0.000000 

# What penalization level leads to good prediction performance?
set.seed(123)
cv_results <- elnet_cv(x, freeny$y, alpha = c(0.5, 0.75),
                       cv_repl = 10, cv_k = 4,
                       cv_measure = "tau")
plot(cv_results, se_mult = 1.5)

plot(cv_results, se_mult = 1.5, what = "coef.path")



# Extract the coefficients at the penalization level with
# smallest prediction error ...
summary(cv_results)
#> EN fit with prediction performance estimated by replications of 4-fold 
#> cross-validation.
#> 
#> 4 out of 4 predictors have non-zero coefficients:
#> 
#>               Estimate
#> (Intercept) -9.6491805
#> X1           0.1899399
#> X2          -0.6858733
#> X3           0.7075924
#> X4           1.2247539
#> ---
#> 
#> Hyper-parameters: lambda=0.003787891, alpha=0.5
coef(cv_results)
#>           (Intercept) lag.quarterly.revenue           price.index 
#>            -9.6491805             0.1899399            -0.6858733 
#>          income.level      market.potential 
#>             0.7075924             1.2247539 
# ... or at the penalization level with prediction error
# statistically indistinguishable from the minimum.
summary(cv_results, lambda = "1.5-se")
#> EN fit with prediction performance estimated by replications of 4-fold 
#> cross-validation.
#> 
#> 4 out of 4 predictors have non-zero coefficients:
#> 
#>               Estimate
#> (Intercept) -9.5875726
#> X1           0.2270959
#> X2          -0.6216322
#> X3           0.6519060
#> X4           1.1972788
#> ---
#> 
#> Hyper-parameters: lambda=0.01303176, alpha=0.5
coef(cv_results, lambda = "1.5-se")
#>           (Intercept) lag.quarterly.revenue           price.index 
#>            -9.5875726             0.2270959            -0.6216322 
#>          income.level      market.potential 
#>             0.6519060             1.1972788 
```
