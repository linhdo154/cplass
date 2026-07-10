# cplass

**CPLASS** (Continuous Piecewise-Linear Approximation via Stochastic
Search) is an R package for detecting changes in velocity within
multidimensional time series, with a focus on intracellular transport
trajectories.

Unlike traditional changepoint methods that target changes in *mean*,
CPLASS models a *continuous* piecewise-linear anchor trajectory and
searches the space of changepoint configurations with a tailored
Metropolis-Hastings sampler. A domain-informed speed penalty discourages
biophysically implausible segment speeds, and the Cumulative Speed
Allocation (CSA) statistic summarizes inferred transport behavior
without requiring an arbitrary motility threshold.

Methodological details, simulation studies, and applications to
lysosomal transport, quantum-dot motor assays, and EB1-GFP microtubule
comet tracking are in the companion paper:

> Do, L., Do, D., Cook, K. J., Shen, Y., & McKinley, S. A. (2026).
> *Change-in-velocity detection in multidimensional data.*

## Installation

CPLASS is not yet on CRAN. Install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("linhdo154/cplass")
```

## Quick start

``` r

library(cplass)

# A small built-in example: lysosome trajectories from the cell periphery
data_file <- system.file("extdata", "Real_21_Periphery.csv", package = "cplass")
traj <- read.csv(data_file)
path <- traj[traj$index_path == 3, ]

fit <- CPLASS(path$t, path$x, path$y, iter_max = 5000, burn_in = 500)

fit$segments_inferred          # inferred segments: duration, speed, state
plot_path_inferred(fit, motor = "Lysosome", max_speed = 2, path_index = 3)
```

The penalized criterion can have multiple local maxima (see Figure 4 of
the paper), so for anything going into a reported result, prefer several
independent starts:

``` r

ms <- CPLASS_multistart(path$t, path$x, path$y, n_starts = 5)
ms$best$segments_inferred
ms$score_spread   # small => starts agree on the same optimum
```

Summarize a fitted trajectory with the Cumulative Speed Allocation
instead of an arbitrary motile/stationary speed cutoff:

``` r

speed_mesh <- seq(0, max(ms$best$segments_inferred$speeds), length.out = 100)
csa <- compute_csa(ms$best$segments_inferred, speed_mesh)
plot(csa$s, csa$csa, type = "l", xlab = "speed", ylab = "cumulative time")
```

See
[`vignette("cplass_example")`](https://linhdo154.github.io/cplass/articles/cplass_example.md)
for a complete walkthrough, or the [package
website](https://linhdo154.github.io/CPLASS/) once published.

## What’s in the box

| Function | Purpose |
|----|----|
| [`CPLASS()`](https://linhdo154.github.io/cplass/reference/CPLASS-function.md) | Run the full algorithm on one 2D trajectory |
| [`CPLASS_multistart()`](https://linhdo154.github.io/cplass/reference/CPLASS_multistart.md) | Multiple independent chains, keep the best |
| [`check_convergence()`](https://linhdo154.github.io/cplass/reference/check_convergence.md) / [`plot_convergence()`](https://linhdo154.github.io/cplass/reference/plot_convergence.md) | Diagnose whether the search has stabilized |
| [`compute_csa()`](https://linhdo154.github.io/cplass/reference/compute_csa.md) / [`csa_theoretical()`](https://linhdo154.github.io/cplass/reference/csa_theoretical.md) | Cumulative Speed Allocation summary |
| [`plot_path_inferred()`](https://linhdo154.github.io/cplass/reference/plot_path_inferred.md) | Four-panel dashboard for one fitted trajectory |
| [`piecewise_linear_con()`](https://linhdo154.github.io/cplass/reference/piecewise_linear_con.md) | Continuous piecewise-linear MLE given changepoints |
| [`MHsearch()`](https://linhdo154.github.io/cplass/reference/MHsearch.md) | The Metropolis-Hastings search itself, for advanced use |

Run [`library(help = "cplass")`](https://rdrr.io/r/base/library.html)
for the full function index.

## Citation

If you use CPLASS in published work, please cite:

    Do L, Do D, Cook KJ, Shen Y, McKinley SA (2026). Change-in-velocity
    detection in multidimensional data. 

A BibTeX entry is available via `citation("cplass")` once installed.

## Contributing / issues

Bug reports and feature requests are welcome at
<https://github.com/linhdo154/cplass/issues>.

## License

MIT © Linh Do and contributors. See
[LICENSE](https://linhdo154.github.io/cplass/LICENSE).
