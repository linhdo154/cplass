# Broadcast per-segment states back onto the original observation grid

Broadcast per-segment states back onto the original observation grid

## Usage

``` r
build_j_inferred(path_inf, seg_inf)
```

## Arguments

- path_inf:

  a \`path_inferred\` tibble (must have column \`t\`).

- seg_inf:

  a \`segments_inferred\` tibble (must have columns \`cp_times\`,
  \`states\`).

## Value

integer vector of per-observation state labels, length
\`nrow(path_inf)\`.
