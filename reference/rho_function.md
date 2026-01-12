# List Available Rho Functions

List Available Rho Functions

## Usage

``` r
rho_function(rho, convex_ok = TRUE)
```

## Arguments

- rho:

  the name of the \\\rho\\ function to check for existence.

- convex_ok:

  if convex \\\rho\\ function is acceptable or not.

## Value

if `rho` is missing returns a vector of supported \\\rho\\ function
names, otherwise the internal integer representation of the \\\rho\\
function.

## See also

Other Robustness control options:
[`consistency_const()`](rho-tuning-constants.md),
[`mscale_algorithm_options()`](mscale_algorithm_options.md)
