# Use the ADMM Elastic Net Algorithm

Use the ADMM Elastic Net Algorithm

## Usage

``` r
en_admm_options(max_it = 1000, step_size, acceleration = 1)
```

## Arguments

- max_it:

  maximum number of iterations.

- step_size:

  step size for the algorithm.

- acceleration:

  acceleration factor for linearized ADMM.

## Value

options for the ADMM EN algorithm.

## See also

Other LS-EN algorithm options:
[`en_algorithm_options`](en_algorithm_options.md),
[`en_cd_options()`](en_cd_options.md),
[`en_dal_options()`](en_dal_options.md),
[`en_lars_options()`](en_lars_options.md)
