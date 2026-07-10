# Gaussian log-likelihood from residual sum of squares

Gaussian log-likelihood from residual sum of squares

## Usage

``` r
loglikelihood(n, RSS, sd = NA)
```

## Arguments

- n:

  number of observations.

- RSS:

  residual sum of squares (summed over both coordinates).

- sd:

  optional known noise standard deviation. If `NA` (default), the noise
  variance is profiled out (concentrated log-likelihood), which is the
  branch used everywhere in
  [`CS`](https://linhdo154.github.io/cplass/reference/CS.md).

## Value

numeric log-likelihood value.
