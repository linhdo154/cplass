# Plot Cumulative Speed Allocation curves

Plot Cumulative Speed Allocation curves

## Usage

``` r
plot_csa(csa, max_error = 0.1, legend = TRUE)
```

## Arguments

- csa:

  a tibble with columns \`s\`, \`csa\`, \`label\`, \`subsample\`,
  \`Hz\`.

- max_error:

  unused, kept for interface compatibility.

- legend:

  logical; overlay the "True" and first subsample curves with a legend.
  Default \`TRUE\`.
