# Continuous piecewise-linear fit given a set of changepoints

Given observed times `t` and coordinates `x`, `y`, and a vector of
changepoint indices `this_cp`, compute the maximum likelihood continuous
piecewise-linear anchor path (Section 2.2 of the companion paper).

## Usage

``` r
piecewise_linear_con(t, x, y, this_cp, old_version = TRUE)
```

## Arguments

- t:

  numeric vector of observation times.

- x:

  numeric vector of x-coordinates.

- y:

  numeric vector of y-coordinates.

- this_cp:

  integer vector of changepoint indices.

- old_version:

  logical. If `TRUE` (default), return the flat list format used
  internally by
  [`CS`](https://linhdo154.github.io/cplass/reference/CS.md) and
  [`MHsearch`](https://linhdo154.github.io/cplass/reference/MHsearch.md).
  If `FALSE`, return the \`tibble\`-based format (\`path_inferred\`,
  \`segments_inferred\`) used by the plotting and CSA helpers.

## Value

See Details; matches the original function's return format.

## Details

This is a drop-in replacement for the original research-script function
of the same name: the argument list and the structure of the returned
list are unchanged, so existing scripts that call
`piecewise_linear_con()` do not need to be modified. Internally it now
delegates to an internal fast fitting routine, which fits both
coordinates from a single QR decomposition instead of two.
