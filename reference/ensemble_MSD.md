# Ensemble (group-averaged) MSD

Ensemble (group-averaged) MSD

## Usage

``` r
ensemble_MSD(msd_df)
```

## Arguments

- msd_df:

  a tibble stacked across paths (e.g. via
  `dplyr::bind_rows(msd_list, .id = "id")`), with columns \`t_lag\` and
  \`msd\`.

## Value

a tibble with columns \`t_lag\` and \`e_msd\`.

## Details

NOTE: the original script's \`ensemble_MSD()\` grouped by a column named
\`t\`, but was always called on plain lists of \`pathwise_MSD()\`
outputs (which use \`t_lag\`, not \`t\`/\`id\`) rather than on a
combined data frame, so \`plot_msd()\`'s calls to
\`ensemble_MSD(msd_list)\` would have errored as originally written.
This version is made internally consistent (grouping by \`t_lag\` on a
properly stacked data frame), but the exact intended grouping/labeling
should be confirmed.
