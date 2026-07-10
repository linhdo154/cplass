# Fill gaps in an unevenly-sampled 2D trajectory

Regularizes a trajectory onto a uniform time grid, linearly
interpolating missing observations (with added noise matched to the
observed increment scale) and replacing implausible jumps.

## Usage

``` r
fill_missing_data(t, x, y, track_id)
```

## Arguments

- t, x, y:

  observed times and coordinates (may have gaps/dropped frames).

- track_id:

  an identifier carried through to the output.

## Value

a list with the regularized \`t\`, \`x\`, \`y\`, indices that were
\`filled\` vs. \`replaced\` for each coordinate, and \`track_id\`.
