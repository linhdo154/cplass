# Alternate 3-row layout for a single segmented trajectory

Alternate 3-row layout for a single segmented trajectory

## Usage

``` r
plot_path_inferred_second_version(
  path_info,
  xy_width = NA,
  t_lim = NA,
  motor,
  max_speed,
  state_shaded = TRUE,
  title_ind = TRUE,
  show_time_changes = FALSE,
  path_index = NULL
)
```

## Arguments

- path_info:

  CPLASS output for one path (\`path_inferred\`, \`segments_inferred\`).

- xy_width:

  width of the x-y panel; computed automatically if \`NA\`.

- t_lim:

  time-axis limits; computed automatically if \`NA\`.

- motor:

  label used in the plot title.

- max_speed:

  upper limit for the speed axis in the duration/speed panel.

- state_shaded:

  logical; shade motile/stationary segments. Default \`TRUE\`.

- title_ind:

  logical; include a path index in the title. Default \`TRUE\`.

- show_time_changes:

  logical; draw vertical lines at changepoints instead of state shading.
  Default \`FALSE\`.

- path_index:

  optional explicit path index/label for the title (see Details in
  package overview about the original's implicit \`i\` lookup).
