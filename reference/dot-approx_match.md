# Approximate Value Matching

Approximate Value Matching

## Usage

``` r
.approx_match(x, table, eps)
```

## Arguments

- x, table:

  see [base::match](https://rdrr.io/r/base/match.html) for details.

- eps:

  numerical tolerance for matching.

## Value

a vector the same length as `x` with integers giving the position in
`table` of the first match if there is a match, or `NA_integer_`
otherwise.
