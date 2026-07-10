# Plot a CPLASS chain's convergence trace

Plot a CPLASS chain's convergence trace

## Usage

``` r
plot_convergence(conv)
```

## Arguments

- conv:

  the output of
  [`check_convergence`](https://linhdo154.github.io/cplass/reference/check_convergence.md).

## Value

a ggplot object: the accepted-state criterion score per iteration
(light) with the running maximum overlaid (dark), and a dashed vertical
line marking the last improvement.
