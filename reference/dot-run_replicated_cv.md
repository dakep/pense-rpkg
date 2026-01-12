# Run replicated K-fold CV with random splits

Run replicated K-fold CV with random splits

## Usage

``` r
.run_replicated_cv(
  std_data,
  cv_k,
  cv_repl,
  cv_est_fun,
  metric,
  par_cluster = NULL,
  handler_args = list()
)
```

## Arguments

- std_data:

  standardized full data set (standardized by `.standardize_data`)

- cv_k:

  number of folds per CV split

- cv_repl:

  number of CV replications.

- cv_est_fun:

  function taking the standardized training set and the indices of the
  left-out observations and returns a list of estimates. The function
  always needs to return the same number of estimates!

- metric:

  function taking a vector of prediction errors and returning the scale
  of the prediction error.

- par_cluster:

  parallel cluster to parallelize computations.

- handler_args:

  additional arguments to the handler function.
