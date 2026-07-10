# Theoretical CSA under a two-state (stationary/motile) kinetic model

Theoretical CSA under a two-state (stationary/motile) kinetic model

## Usage

``` r
csa_theoretical(theta, s_mesh)
```

## Arguments

- theta:

  a list with elements \`speed_alpha\`, \`speed_beta\`, \`dist_avg\`,
  \`p_MM\`, \`dur_stationary_avg\`, \`DEPENDENT_SPEED_DUR\`.

- s_mesh:

  numeric vector of speeds.

## Value

numeric vector, the theoretical CSA evaluated at \`s_mesh\`.
