# cplass: Continuous Piecewise-Linear Approximation via Stochastic Search

CPLASS detects changes in velocity ("change-in-slope") in
multidimensional time series such as intracellular particle
trajectories. It models observed positions as Gaussian fluctuations
around a continuous piecewise-linear anchor path, and searches the space
of changepoint configurations using a Metropolis-Hastings sampler with
proposal moves tailored to the change-in-velocity problem (see Section
2.4 of the companion paper).

## Details

The main entry points are:

- [`CPLASS`](https://linhdo154.github.io/cplass/reference/CPLASS-function.md):
  run the full algorithm on a single 2D trajectory (t, x, y).

- [`compute_csa`](https://linhdo154.github.io/cplass/reference/compute_csa.md),
  [`summarize_segments_inferred`](https://linhdo154.github.io/cplass/reference/summarize_segments_inferred.md):
  summarize inferred segmentations with the Cumulative Speed Allocation
  statistic.

- [`plot_path_inferred`](https://linhdo154.github.io/cplass/reference/plot_path_inferred.md):
  visualize a single segmented trajectory.

## Performance notes

The internal fitting routine (`piecewise_linear_con`) and the MCMC
driver (`MHsearch`) were rewritten for speed relative to the original
research-script implementation, while keeping identical function
signatures, return structures, and statistical behavior. See
`inst/benchmarks/` (in the GitHub source repo) for benchmark details and
a full list of changes.

## See also

Useful links:

- <https://github.com/linhdo154/cplass>

- Report bugs at <https://github.com/linhdo154/cplass/issues>

## Author

**Maintainer**: Linh Do <thuylinh.do@duke.edu>

Authors:

- Scott A. McKinley
