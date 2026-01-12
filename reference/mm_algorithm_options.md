# MM-Algorithm to Compute Penalized Elastic Net S- and M-Estimates

Additional options for the MM algorithm to compute EN S- and
M-estimates.

## Usage

``` r
mm_algorithm_options(
  max_it = 500,
  tightening = c("adaptive", "exponential", "none"),
  tightening_steps = 2,
  en_algorithm_opts
)
```

## Arguments

- max_it:

  maximum number of iterations.

- tightening:

  how to make inner iterations more precise as the algorithm approaches
  a local minimum.

- tightening_steps:

  for *adaptive* tightening strategy, how often to tighten until the
  desired tolerance is attained.

- en_algorithm_opts:

  options for the inner LS-EN algorithm. See
  [en_algorithm_options](en_algorithm_options.md) for details.

## Value

options for the MM algorithm.

## See also

cd_algorithm_options for a direct optimization of the non-convex PENSE
loss.

Other Robust EN algorithms:
[`cd_algorithm_options()`](cd_algorithm_options.md)
