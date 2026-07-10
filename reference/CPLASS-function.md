# Run CPLASS on a single 2D trajectory

Fits the CPLASS model (Section 2 of the companion paper) to a single
particle trajectory given by observed times `t` and coordinates `x`,
`y`: runs the Metropolis-Hastings search
([`MHsearch`](https://linhdo154.github.io/cplass/reference/MHsearch.md))
over changepoint configurations, selects the configuration with the
highest criterion score, and returns the corresponding continuous
piecewise-linear segmentation.

## Usage

``` r
CPLASS(
  t,
  x,
  y,
  lambda_r = 1/30,
  iter_max = 5000,
  burn_in = 500,
  s_cap = 1,
  gamma = 1.01,
  speed_pen = TRUE,
  Diagnostic = FALSE,
  sd = NA,
  pen = "ssic",
  show_progress = TRUE,
  patience = NULL
)
```

## Arguments

- t:

  numeric vector of observation times.

- x:

  numeric vector of x-coordinates.

- y:

  numeric vector of y-coordinates.

- lambda_r:

  rate parameter for the Type-1 (independent) changepoint proposal.
  Default `1/30`.

- iter_max:

  number of MCMC iterations. Default 5000.

- burn_in:

  number of initial iterations discarded. Default 500.

- s_cap:

  maximum segment speed with no penalty. Default 1.

- gamma:

  exponent for the sSIC penalty. Default 1.01 (see Section 3.1.1 of the
  companion paper for the rationale).

- speed_pen:

  logical; whether to activate the speed penalty. Default `TRUE`.

- Diagnostic:

  logical; if `TRUE`, also return the MCMC diagnostic trace
  (\`info_table\`, \`update_info\`). Default `FALSE`.

- sd:

  optional known noise sd.

- pen:

  one of `"ssic"`, `"aicc"`, `"hybrid"`.

- show_progress:

  logical; display a text progress bar. Default `TRUE`.

- patience:

  if not `NULL` (default), stop the search early once the criterion
  score hasn't improved for this many consecutive iterations, rather
  than always running the full \`iter_max\`. See
  [`MHsearch`](https://linhdo154.github.io/cplass/reference/MHsearch.md)'s
  "Early stopping" section for details and caveats (in particular: this
  doesn't protect against local maxima — pair it with
  [`CPLASS_multistart`](https://linhdo154.github.io/cplass/reference/CPLASS_multistart.md)
  if that's a concern).

## Value

A list with \`segments_inferred\` and \`path_inferred\` tibbles (see
[`piecewise_linear_con`](https://linhdo154.github.io/cplass/reference/piecewise_linear_con.md),
`old_version = FALSE`), plus \`best_score\`: the criterion value
\\\Phi(r)\\ of the selected configuration, useful for comparing chains
(see
[`CPLASS_multistart`](https://linhdo154.github.io/cplass/reference/CPLASS_multistart.md),
[`check_convergence`](https://linhdo154.github.io/cplass/reference/check_convergence.md)),
and \`stopped_early\` / \`iterations_used\` (only meaningful if
\`patience\` was set). If `Diagnostic = TRUE`, also includes
\`info_table\` and \`update_info\`. Returns `pl = NULL` (wrapped the
same way) if the sampler fails to converge on a valid fit after retries.

## Details

Same signature and return value as the original research-script
`CPLASS()`. See
[`MHsearch`](https://linhdo154.github.io/cplass/reference/MHsearch.md)
for details on what changed internally for speed.

## Examples

``` r
# \donttest{
data_file <- system.file("extdata", "Real_21_Periphery.csv", package = "cplass")
traj <- read.csv(data_file)
path <- traj[traj$index_path == 3, ]

fit <- CPLASS(path$t, path$x, path$y, iter_max = 2000, burn_in = 200,
               show_progress = FALSE)
fit$segments_inferred
#> # A tibble: 5 × 7
#>   cp_times durations states speeds       vx      vy angles
#>      <dbl>     <dbl>  <int>  <dbl>    <dbl>   <dbl>  <dbl>
#> 1     10.5    10.4        0 0.0286 -0.00458  0.0283  1.73 
#> 2     16.2     5.65       0 0.0529  0.00748 -0.0523 -1.43 
#> 3     23.2     7          0 0.0136 -0.00638  0.0120  2.06 
#> 4     23.6     0.450      1 0.440   0.157   -0.411  -1.21 
#> 5     30.0     6.35       0 0.0164  0.0105  -0.0127 -0.881
# }
```
