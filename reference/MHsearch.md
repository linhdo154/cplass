# Metropolis-Hastings search over the changepoint space

Runs the tailored MH sampler described in Section 2.4 of the companion
paper for a single 2D trajectory. This has the same signature and return
structure as the original research-script `MHsearch()`, but is
substantially faster internally. See "Performance" below.

## Usage

``` r
MHsearch(
  t,
  x,
  y,
  dt,
  lambda_r = 1/30,
  iter_max = 5000,
  burn_in = 500,
  s_cap = 1,
  gamma = 1.01,
  speed_control = 0,
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

- dt:

  observation time step.

- lambda_r:

  rate parameter for the Type-1 (independent) changepoint proposal.

- iter_max:

  number of MCMC iterations (an upper bound if \`patience\` is set; see
  "Early stopping").

- burn_in:

  number of initial iterations discarded from the trace.

- s_cap:

  maximum segment speed with no penalty.

- gamma:

  exponent for the sSIC penalty.

- speed_control:

  1 to activate the speed penalty, 0 to deactivate.

- sd:

  optional known noise sd (see
  [`loglikelihood`](https://linhdo154.github.io/cplass/reference/loglikelihood.md)).

- pen:

  one of `"ssic"`, `"aicc"`, `"hybrid"`.

- show_progress:

  logical; display a text progress bar (default `TRUE`, matching the
  original behavior).

- patience:

  if not `NULL` (default), stop early once the post-burn-in
  running-maximum criterion score hasn't improved for this many
  consecutive iterations. See "Early stopping".

## Value

A list with \`cps_list\` (post-burn-in trace), \`info_table\`
(post-burn-in diagnostics), \`update_info\` (acceptance rate), and
\`stopped_early\` / \`iterations_used\` (whether/when patience
triggered; \`iterations_used\` equals \`iter_max\` if it never did).

## Performance

The original implementation recomputed the full piecewise-linear fit
([`piecewise_linear_con`](https://linhdo154.github.io/cplass/reference/piecewise_linear_con.md))
up to three times per iteration: once for the current state, once for
the proposal, and once more after the accept/reject step to populate the
trace. Since the fit is a pure function of `(t, x, y, changepoints)`,
the current state's fit never needs to be recomputed (it's identical to
the state fit at the end of the previous iteration), and the
post-decision refit is always identical to either the already-computed
current or proposal fit. This version fits the proposal once per
iteration and carries the accepted fit forward, cutting the number of
\`piecewise_linear_con()\`/\`CS()\` calls from ~3 to 1 per iteration.
The diagnostic trace (\`info_table\`) is also built from preallocated
vectors instead of row-by-row \`bind_rows()\`, which previously made the
trace step scale quadratically in \`iter_max\`.

Statistically, this produces the same target distribution and the same
accept/reject decisions as the original for a given stream of random
draws in the same order. Because the birth-vector proposal now uses
[`stats::rbinom()`](https://rdrr.io/r/stats/Binomial.html) instead of
`Rlab::rbern()`, results will not be bit-identical to the original code
for the same [`set.seed()`](https://rdrr.io/r/base/Random.html), but are
statistically equivalent (both are Bernoulli draws).

## Early stopping

If `patience` is set, `iter_max` becomes a safety ceiling rather than a
fixed target: once the running-maximum criterion score (post-burn-in)
hasn't improved for `patience` consecutive iterations, the chain stops
early. This is the same "no improvement" idea as
[`check_convergence`](https://linhdo154.github.io/cplass/reference/check_convergence.md),
applied live inside the loop so you don't pay for iterations a path
doesn't need. It does not protect against local maxima in the criterion
surface (see
[`CPLASS_multistart`](https://linhdo154.github.io/cplass/reference/CPLASS_multistart.md)
for that) — a chain can plateau on a local max just as easily as a
global one.
