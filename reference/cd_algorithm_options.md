# Coordinate Descent (CD) Algorithm to Compute Penalized Elastic Net S-estimates

Set options for the CD algorithm to compute adaptive EN S-estimates.

## Usage

``` r
cd_algorithm_options(
  max_it = 1000,
  reset_it = 8,
  linesearch_steps = 4,
  linesearch_mult = 0.5
)
```

## Arguments

- max_it:

  maximum number of iterations.

- reset_it:

  number of iterations after which the residuals are re-computed from
  scratch, to prevent numerical drifts from incremental updates.

- linesearch_steps:

  maximum number of steps used for line search.

- linesearch_mult:

  multiplier to adjust the step size in the line search.

## Value

options for the CD algorithm to compute (adaptive) PENSE estimates.

## See also

mm_algorithm_options to optimize the non-convex PENSE objective function
via a sequence of convex problems.

Other Robust EN algorithms:
[`mm_algorithm_options()`](mm_algorithm_options.md)
