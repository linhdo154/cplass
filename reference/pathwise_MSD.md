# Pathwise mean squared displacement (MSD)

Pathwise mean squared displacement (MSD)

## Usage

``` r
pathwise_MSD(path, use_pct = 0.2)
```

## Arguments

- path:

  a path object with \`path\$path\$x\`, \`path\$path\$y\`,
  \`path\$path\$t\`.

- use_pct:

  fraction of the trajectory length to use as the maximum lag. Default
  0.2.

## Value

a tibble with columns \`t_lag\` and \`msd\`.
