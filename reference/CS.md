# Criterion function for a changepoint configuration

Computes \\\Phi(r) = 2 \hat L_n - \mathrm{pen}(r)\\ (Definition 1/2 of
the companion paper): the penalized log-likelihood of the continuous
piecewise-linear fit associated with a changepoint vector `r`, using
either the strengthened SIC penalty, AICc, or the max of the two
("hybrid"), plus an optional speed penalty.

## Usage

``` r
CS(t, x, y, r, p1, s_cap, with_speed, sd = NA, pen = "ssic", gamma)
```

## Arguments

- t, x, y:

  trajectory data.

- r:

  changepoint indicator vector (0/1), length `length(t) - 2`.

- p1:

  penalty coefficient \\(\log n)^\gamma\\ for the linear (sSIC) term;
  only used when `pen == "ssic"` indirectly via the caller's precomputed
  value (kept for interface compatibility with the original script;
  `gamma` is what actually parametrizes sSIC here).

- s_cap:

  maximum segment speed with no penalty.

- with_speed:

  1 to activate the speed penalty, 0 to deactivate.

- sd:

  optional known noise sd, passed to
  [`loglikelihood`](https://linhdo154.github.io/cplass/reference/loglikelihood.md).

- pen:

  one of `"ssic"`, `"aicc"`, `"hybrid"`.

- gamma:

  exponent for the sSIC penalty.

## Value

A list with \`s\` (criterion value), \`llh\`, \`p\` (penalty), \`pv\`
(speed-penalty component), \`logical\`.
