# Side-by-side comparison of true vs. inferred segmentation

Side-by-side comparison of true vs. inferred segmentation

## Usage

``` r
plot_path_actual_and_inferred(
  path_info,
  xy_width = NA,
  t_lim = NA,
  state_shaded = TRUE,
  path_index = NULL
)
```

## Arguments

- path_info:

  a path object with both ground-truth (\`path\`, \`segments\`) and
  inferred (\`path_inferred\`, \`segments_inferred\`) fields.

- xy_width:

  width of the x-y panel; computed automatically if \`NA\`.

- t_lim:

  time-axis limits; computed automatically if \`NA\`.

- state_shaded:

  logical; shade motile/stationary segments. Default \`TRUE\`.

- path_index:

  optional explicit path index/label for the title (see Details in
  package overview about the original's implicit \`i\` lookup).
