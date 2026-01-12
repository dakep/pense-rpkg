# Run replicated K-fold CV with random splits, matching the global estimates to the CV estimates by Kendall's tau-b computed on the robustness weights.

Run replicated K-fold CV with random splits, matching the global
estimates to the CV estimates by Kendall's tau-b computed on the
robustness weights.

## Usage

``` r
.run_replicated_cv_ris(
  std_data,
  cv_k,
  cv_repl,
  cv_est_fun,
  global_ests,
  min_similarity = 0,
  par_cluster = NULL,
  rho_opts,
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

- global_ests:

  estimates computed on all observations.

- min_similarity:

  minimum (average) similarity for CV solutions to be considered
  (between 0 and 1). If no CV solution satisfies this lower bound, the
  best CV solution will be used regardless of similarity.

- par_cluster:

  parallel cluster to parallelize computations.

- rho_opts:

  rho function options.

- handler_args:

  additional arguments to the handler function.
