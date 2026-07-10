# Convert a changepoint index vector to a 0/1 indicator vector

Convert a changepoint index vector to a 0/1 indicator vector

## Usage

``` r
cp_to_r(cp, n)
```

## Arguments

- cp:

  integer vector of changepoint indices.

- n:

  number of observations.

## Examples

``` r
cp_to_r(cp = c(5, 12), n = 20)
#>  [1] 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0
```
