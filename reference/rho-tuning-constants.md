# Get the Constant for Consistency for the M-Scale and for Efficiency for the M-estimate of Location

Returns the tuning constants required to achieve the desired breakdown
point or efficiency under the Normal model.

## Usage

``` r
consistency_const(delta, rho, eps = sqrt(.Machine$double.eps))

efficiency_const(eff, rho, eps = sqrt(.Machine$double.eps))
```

## Arguments

- delta:

  desired breakdown point (between 0 and 0.5)

- rho:

  the name of the chosen \\\rho\\ function. See
  [`rho_function()`](rho_function.md) for a list of supported functions.

- eps:

  numerical tolerance level for equality comparisons

- eff:

  desired asymptotic efficiency (between 0.1 and 0.99).

## Value

consistency constant

## See also

Other Robustness control options:
[`mscale_algorithm_options()`](mscale_algorithm_options.md),
[`rho_function()`](rho_function.md)

Other Robustness control options:
[`mscale_algorithm_options()`](mscale_algorithm_options.md),
[`rho_function()`](rho_function.md)
