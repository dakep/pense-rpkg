# Options for the M-scale Estimation Algorithm

Options for the M-scale Estimation Algorithm

## Usage

``` r
mscale_algorithm_options(rho = "bisquare", max_it = 200, eps = 1e-08)
```

## Arguments

- rho:

  the \\\rho\\ function to use. See [`rho_function()`](rho_function.md)
  for possible values.

- max_it:

  maximum number of iterations.

- eps:

  numerical tolerance to check for convergence.

## Value

options for the M-scale estimation algorithm.

## See also

Other Robustness control options:
[`consistency_const()`](rho-tuning-constants.md),
[`rho_function()`](rho_function.md)
