# Total duration over which true and inferred segmentations disagree on state

Sweeps two segmentations (true and CPLASS-inferred) simultaneously with
a two-pointer merge over their changepoints and sums the overlap
duration where the two disagree on motile/stationary state.

## Usage

``` r
duration_of_difference(path_info)
```

## Arguments

- path_info:

  a path object with \`path\`, \`segments\`, \`segments_inferred\`.

## Value

numeric, total duration of disagreement (same units as \`t\`).
