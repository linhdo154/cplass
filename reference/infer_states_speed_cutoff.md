# Classify inferred segments as motile/stationary by a speed cutoff

Classify inferred segments as motile/stationary by a speed cutoff

## Usage

``` r
infer_states_speed_cutoff(seg_speeds, cutoff = 0.1)
```

## Arguments

- seg_speeds:

  numeric vector of inferred segment speeds.

- cutoff:

  speed threshold; segments faster than this are labeled 1 (motile),
  others 0 (stationary). Default 0.1.

## Value

integer vector of 0/1 labels, same length as \`seg_speeds\`.
