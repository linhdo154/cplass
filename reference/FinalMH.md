# Collect a Metropolis-Hastings trace into flat vectors/lists

Drop-in, vectorized replacement for the original `FinalMH()`. The
original also computed a \`loglikelihood\` vector that was never
included in its return value (dead code); this version omits that
computation.

## Usage

``` r
FinalMH(bcp_new)
```

## Arguments

- bcp_new:

  the \`cps_list\` element of
  [`MHsearch`](https://linhdo154.github.io/cplass/reference/MHsearch.md)'s
  output.

## Value

A list of parallel vectors/lists extracted from the trace: \`r\`, \`u\`,
\`v\`, \`this_cp\`, \`seg_vel\`, \`seg_time\`, \`eta\`,
\`criterion_score\`, \`RSS\`, \`num_cp\`, \`r_prop\`.
