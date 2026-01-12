# Use Coordinate Descent to Solve Elastic Net Problems

Use Coordinate Descent to Solve Elastic Net Problems

## Usage

``` r
en_cd_options(max_it = 1000, reset_it = 8)
```

## Arguments

- max_it:

  maximum number of iterations.

- reset_it:

  number of iterations after which the residuals are re-computed from
  scratch, to prevent numerical drifts from incremental updates.

## See also

Other LS-EN algorithm options:
[`en_admm_options()`](en_admm_options.md),
[`en_algorithm_options`](en_algorithm_options.md),
[`en_dal_options()`](en_dal_options.md),
[`en_lars_options()`](en_lars_options.md)
