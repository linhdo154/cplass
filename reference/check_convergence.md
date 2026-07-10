# Check whether a single CPLASS chain has stopped improving

CPLASS's Metropolis-Hastings search does not converge to a stationary
distribution in the usual Bayesian-inference sense — it is a stochastic
search for the changepoint configuration that maximizes the criterion
\\\Phi(r)\\. The relevant question for a single chain is therefore not
"has it mixed?" but "has it stopped finding better configurations?".

## Usage

``` r
check_convergence(
  t,
  x,
  y,
  lambda_r = 1/30,
  iter_max = 5000,
  burn_in = 500,
  s_cap = 1,
  gamma = 1.01,
  speed_pen = TRUE,
  sd = NA,
  pen = "ssic",
  window_frac = 0.2,
  show_progress = FALSE
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

- sd:

  optional known noise sd.

- pen:

  one of `"ssic"`, `"aicc"`, `"hybrid"`.

- window_frac:

  fraction of post-burn-in iterations used as the "no improvement"
  window. Default 0.2 (last 20%).

- show_progress:

  logical; display a text progress bar. Default `TRUE`.

## Value

A list with:

- \`converged\`: logical.

- \`best_score\`: the best criterion score found.

- \`last_improvement_iter\`: post-burn-in iteration index at which the
  running maximum last improved.

- \`iterations_since_improvement\`, \`window_size\`: the check's inputs,
  for transparency.

- \`trace\`: a tibble with \`iteration\`, \`score\`, \`running_max\`,
  for plotting (see
  [`plot_convergence`](https://linhdo154.github.io/cplass/reference/plot_convergence.md)).

## Details

This runs
[`CPLASS`](https://linhdo154.github.io/cplass/reference/CPLASS-function.md)
once with `Diagnostic = TRUE`, reconstructs the trace of the
\*accepted\* state's criterion score at every iteration (from
\`info_table\`'s \`score_cur\`/\`score_new\`/\`decision\` columns), and
checks how long it has been since the running maximum last improved. If
nothing better has been found in the last `window_frac` fraction of
post-burn-in iterations, the chain is flagged as (locally) converged.

Important: this only tells you the \*chain\* has settled, not that it
found the \*global\* optimum — the criterion surface can have multiple
competing local maxima (see Figure 4 of the companion paper), which is
exactly what
[`CPLASS_multistart`](https://linhdo154.github.io/cplass/reference/CPLASS_multistart.md)
is for. Use both together: this to size \`iter_max\` correctly for a
given trajectory length, and multistart to guard against local-maximum
disagreement between runs.

## Examples

``` r
# \donttest{
data_file <- system.file("extdata", "Real_21_Periphery.csv", package = "cplass")
traj <- read.csv(data_file)
path <- traj[traj$index_path == 3, ]

conv <- check_convergence(path$t, path$x, path$y, iter_max = 2000, burn_in = 200)
conv$converged
#> [1] FALSE
# }
```
