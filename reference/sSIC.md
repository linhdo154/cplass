# Strengthened Schwarz Information Criterion penalty

\\(\log n)^\gamma \rho\\, as in Definition 2 of the companion paper.

## Usage

``` r
sSIC(ncp, n, gamma)
```

## Arguments

- ncp:

  number of changepoints.

- n:

  number of observations.

- gamma:

  exponent, `gamma > 1` recommended (default in
  [`CPLASS`](https://linhdo154.github.io/cplass/reference/CPLASS-function.md)
  is 1.01).

## Examples

``` r
sSIC(ncp = 2, n = 200, gamma = 1.01)
#> [1] 59.26142
```
