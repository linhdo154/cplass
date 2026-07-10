# Proposal functions for the CPLASS Metropolis-Hastings sampler

These implement the four proposal types described in Section 2.4 of the
companion paper: (1) an independent changepoint draw (\`q_new\`), (2)
birth/death of a single changepoint (\`q_bd\`), (3) birth/death of a
segment / paired changepoints (\`q_bd2\`, via \`q_as\`/\`q_ds\`), and
(4) a location shift of a single changepoint (\`q_shift\`).
\`proposal_function()\` combines them with mixture weights \`u1, u2,
u3\`; \`pproposal()\` evaluates the corresponding proposal density for
the Metropolis-Hastings ratio.

## Usage

``` r
log_q_new_pmf(lambda_r, dt, N, r)

q_new(lambda_r, dt, N)

log_q_bd_pmf(r, status)

q_bd(r_cur)

q_shift(r_cur)

q_as_pmf(r)

q_as(r_cur)

q_ds_pmf(r)

q_ds(r_cur)

q_bd2(r_cur)

log_q_bd2_pmf(r, status)

proposal_function(u, r_cur, N, lambda_r, dt)

pproposal(u, r, N, lambda_r, r_given, status, dt)
```

## Arguments

- lambda_r:

  rate parameter for the Type-1 (independent) proposal.

- dt:

  observation time step.

- N:

  number of observations.

- status:

  proposal-type-specific status flag returned alongside \`r_prop\` by
  \`q_bd()\`/\`q_bd2()\`.

- r_cur, r:

  current changepoint indicator vector.

- u:

  draw from Uniform(0,1) selecting which proposal type to use.

- r_given:

  the changepoint vector whose proposal density is being evaluated
  (either the proposed or current vector, depending on direction of the
  MH ratio).
