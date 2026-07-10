# Convert a 0/1 changepoint indicator vector to changepoint indices

Convert a 0/1 changepoint indicator vector to changepoint indices

## Usage

``` r
r_to_cp(r)
```

## Arguments

- r:

  indicator vector.

## Examples

``` r
r_to_cp(cp_to_r(cp = c(5, 12), n = 20))
#> [1]  5 12
```
