# Run CPLASS from multiple independent random starts and keep the best

Runs
[`CPLASS`](https://linhdo154.github.io/cplass/reference/CPLASS-function.md)
`n_starts` times with different seeds and returns the run with the
highest criterion score (`best_score`), along with a summary of how much
the runs agreed — the multi-chain analogue of
[`check_convergence`](https://linhdo154.github.io/cplass/reference/check_convergence.md)'s
single-chain check. Large disagreement in `all_scores` or
`all_n_segments` indicates the criterion surface has multiple
comparably-good local maxima for this trajectory (see Figure 4 of the
companion paper) — in that case, the \*Cumulative Speed Allocation\*
summary
([`compute_csa`](https://linhdo154.github.io/cplass/reference/compute_csa.md))
is a more robust downstream quantity than any single run's segment
count.

## Usage

``` r
CPLASS_multistart(
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
  n_starts = 5,
  patience = "auto",
  seeds = NULL,
  n_cores = 1,
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

- n_starts:

  number of independent chains to run. Default 5.

- patience:

  early-stopping patience passed to each chain's
  [`MHsearch`](https://linhdo154.github.io/cplass/reference/MHsearch.md)
  call. Default `"auto"`: uses `ceiling(0.3 * iter_max)`. Set to `NULL`
  to disable early stopping (every chain runs the full \`iter_max\`), or
  supply a number to use that patience directly. See "Patience +
  multistart, together".

- seeds:

  optional integer vector of length \`n_starts\`; if \`NULL\` (default),
  seeds are drawn from the current RNG state, and are returned in the
  output so the winning run can be reproduced exactly with
  `set.seed(result$best_seed)`.

- n_cores:

  number of cores to use via
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (not available on Windows; falls back to sequential with a message).
  Default 1 (sequential).

- show_progress:

  logical; display a text progress bar. Default `TRUE`.

## Value

A list with:

- \`best\`: the winning
  [`CPLASS`](https://linhdo154.github.io/cplass/reference/CPLASS-function.md)
  output (highest \`best_score\`).

- \`best_seed\`: the seed that produced it — pass to \`set.seed()\`
  before re-running \`CPLASS()\` to reproduce it exactly.

- \`all_scores\`, \`all_n_segments\`, \`seeds\`: per-chain results, for
  inspecting agreement across starts.

- \`score_spread\`: \`max(all_scores) - min(all_scores)\`; small values
  indicate the starts agree on (approximately) the same optimum.

- \`patience_used\`: the actual patience value applied to every chain
  (resolved from \`"auto"\` if applicable), for transparency.

## Patience + multistart, together

By default (`patience = "auto"`), each of the `n_starts` chains uses
[`MHsearch`](https://linhdo154.github.io/cplass/reference/MHsearch.md)'s
patience-based early stopping (ceiling(0.3 \* iter_max) iterations
without improvement), so no single chain pays for iterations it doesn't
need. Multistart's job is to make that safe: because each chain has an
independent random trajectory through the search space, it's very
unlikely that two different chains plateau at exactly the same (possibly
suboptimal) configuration, so taking the max across chains recovers what
any one early-stopped chain might have missed. Empirically (see the
package's benchmark scripts), \`n_starts = 5\` with auto patience
matches the best score of 5 full-length (no-patience) runs almost
exactly, at a fraction of the time. Set `patience = NULL` to disable
early stopping entirely and run every chain to the full \`iter_max\`
(slower, and in our tests did not find better results than the default —
but available if you want to verify that yourself on your own data).

## Examples

``` r
# \donttest{
data_file <- system.file("extdata", "Real_21_Periphery.csv", package = "cplass")
traj <- read.csv(data_file)
path <- traj[traj$index_path == 3, ]

ms <- CPLASS_multistart(path$t, path$x, path$y, iter_max = 2000, burn_in = 200,
                         n_starts = 3, show_progress = FALSE)
ms$all_scores
#> [1] 2143.638 2146.754 2176.594
ms$score_spread
#> [1] 32.95607
# }
```
