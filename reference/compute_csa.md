# Cumulative Speed Allocation (CSA)

For every speed in `speed_mesh`, computes the duration-weighted
proportion of time spent at speeds less than or equal to that value
(Cook et al. 2025; Section 3.1.2 of the companion paper).

## Usage

``` r
compute_csa(segments_summary, speed_mesh)
```

## Arguments

- segments_summary:

  a data frame/tibble with columns \`durations\` and \`speeds\`.

- speed_mesh:

  numeric vector of speed values at which to evaluate the CSA.

## Value

a tibble with columns \`s\` (the mesh) and \`csa\` (cumulative
proportion of time).

## Examples

``` r
segs <- data.frame(durations = c(2, 1, 3), speeds = c(0.01, 0.5, 0.05))
compute_csa(segs, speed_mesh = seq(0, 0.6, by = 0.1))
#> # A tibble: 7 × 2
#>       s   csa
#>   <dbl> <dbl>
#> 1   0   0    
#> 2   0.1 0.833
#> 3   0.2 0.833
#> 4   0.3 0.833
#> 5   0.4 0.833
#> 6   0.5 0.833
#> 7   0.6 1    
```
