# Just the x-y panel of a segmented trajectory

Just the x-y panel of a segmented trajectory

## Usage

``` r
plot_path_inferred_xy(
  path_info,
  xy_width = NA,
  t_lim = NA,
  motor,
  max_speed,
  state_shaded = TRUE
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
