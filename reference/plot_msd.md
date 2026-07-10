# Plot pathwise and ensemble MSD curves for a group, split by activity

Plot pathwise and ensemble MSD curves for a group, split by activity

## Usage

``` r
plot_msd(msd_list, is_active, group_name = "")
```

## Arguments

- msd_list:

  a list of per-path MSD tibbles (\`pathwise_MSD()\` output).

- is_active:

  logical vector, same length as \`msd_list\`.

- group_name:

  label used in the plot title.

## Value

invisibly, the overall ensemble MSD tibble.
