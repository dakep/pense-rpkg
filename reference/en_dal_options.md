# Use the DAL Elastic Net Algorithm

Use the DAL Elastic Net Algorithm

## Usage

``` r
en_dal_options(
  max_it = 100,
  max_inner_it = 100,
  eta_multiplier = 2,
  eta_start_conservative = 0.01,
  eta_start_aggressive = 1,
  lambda_relchange_aggressive = 0.25
)
```

## Arguments

- max_it:

  maximum number of (outer) iterations.

- max_inner_it:

  maximum number of (inner) iterations in each outer iteration.

- eta_multiplier:

  multiplier for the barrier parameter. In each iteration, the barrier
  must be more restrictive (i.e., the multiplier must be \> 1).

- eta_start_conservative:

  conservative initial barrier parameter. This is used if the previous
  penalty is undefined or too far away.

- eta_start_aggressive:

  aggressive initial barrier parameter. This is used if the previous
  penalty is close.

- lambda_relchange_aggressive:

  how close must the lambda parameter from the previous penalty term be
  to use an aggressive initial barrier parameter (i.e., what constitutes
  "too far").

## Value

options for the DAL EN algorithm.

## See also

Other LS-EN algorithm options:
[`en_admm_options()`](en_admm_options.md),
[`en_algorithm_options`](en_algorithm_options.md),
[`en_cd_options()`](en_cd_options.md),
[`en_lars_options()`](en_lars_options.md)
