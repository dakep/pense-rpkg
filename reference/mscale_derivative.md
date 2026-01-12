# Compute the Gradient and Hessian of the M-Scale Function

Compute the derivative (gradient) or the Hessian of the M-scale function
evaluated at the point `x`.

## Usage

``` r
mscale_derivative(
  x,
  bdp = 0.25,
  order = 1,
  cc,
  opts = mscale_algorithm_options()
)
```

## Arguments

- x:

  numeric values. Missing values are verbosely ignored.

- bdp:

  desired breakdown point (between 0 and 0.5).

- order:

  compute the gradient (`order=1`) or the gradient and the Hessian
  (`order=2`).

- cc:

  cutoff value for the bisquare rho function. By default, chosen to
  yield a consistent estimate for the Normal distribution.

- opts:

  a list of options for the M-scale estimation algorithm, see
  [`mscale_algorithm_options()`](mscale_algorithm_options.md) for
  details.

## Value

a vector of derivatives of the M-scale function, one per element in `x`.
