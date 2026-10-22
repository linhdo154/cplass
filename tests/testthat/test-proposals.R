test_that("case 9: outer move-type probabilities are always 1/4,1/8,1/8,1/2, at all states", {
  # Note: the outer move-type draw is provably state-independent by
  # construction (proposal_function() draws u and picks a type BEFORE
  # doing anything with r_cur), so this is really confirming there's no
  # accidental state-dependence leak, not testing something that could
  # plausibly be close-but-wrong. Absolute-difference checks (not
  # testthat's `tolerance`, whose relative semantics for vectors are not
  # straightforwardly per-element -- see revision notes) with a wide
  # margin relative to the ~0.002-0.003 expected SE at nsim=20000.
  set.seed(101)
  check_weights <- function(r_cur, N, K_max = NULL, nsim = 20000, label = "") {
    counts <- integer(4)
    for (i in seq_len(nsim)) {
      pp <- proposal_function(r_cur, N, lambda_r = 1 / 30, dt = 0.05, K_max = K_max)
      counts[pp$move_type] <- counts[pp$move_type] + 1L
    }
    freq <- counts / nsim
    max_dev <- max(abs(freq - c(0.25, 0.125, 0.125, 0.5)))
    expect_lt(max_dev, 0.02)
  }

  N <- 20
  check_weights(integer(N - 2), N, label = "all-zero state")

  ordinary <- integer(N - 2); ordinary[c(5, 12)] <- 1L
  check_weights(ordinary, N, label = "ordinary interior state")

  near_kmax <- integer(N - 2); near_kmax[1:5] <- 1L # K(r) = 6, near K_max = 6
  check_weights(near_kmax, N, K_max = 6, label = "near K_max")

  saturated_boundary <- integer(N - 2); saturated_boundary[1:(N - 3)] <- 1L # c(r) = n-3, at structural ceiling
  check_weights(saturated_boundary, N, label = "at structural c(r) = n-3 ceiling")
})

test_that("case 10: Type 2 unavailable direction is a self-transition, not redistribution", {
  set.seed(102)
  N <- 20
  all_zero <- integer(N - 2) # death is impossible (c(r) = 0)

  self_loops <- 0L
  births <- 0L
  nsim <- 20000
  for (i in seq_len(nsim)) {
    rp <- .q_bd_core(all_zero, N)$r_prop
    if (identical(rp, all_zero)) self_loops <- self_loops + 1L else births <- births + 1L
  }
  # exactly the "death" half of draws should self-loop (birth should almost
  # always succeed from an all-zero state with plenty of room); confirm the
  # self-loop rate is close to 0.5, not near 0 (which would indicate birth
  # was silently substituted for the blocked death direction) or near 1.
  expect_equal(self_loops / nsim, 0.5, tolerance = 0.03)

  # the log-pmf self-loop mass should independently predict the same ~0.5
  self_mass <- exp(log_q_bd_pmf(all_zero, all_zero, N))
  expect_equal(self_mass, 0.5, tolerance = 1e-8)
})

test_that("case 11: Type 3 insertion never samples an existing segment boundary", {
  N <- 25
  r <- integer(N - 2)
  r[c(5, 15)] <- 1L
  boundaries <- c(1L, r_to_cp(r), N)

  # the closed-form capacity/sampling mechanism never materializes a pair
  # list, so verify the invariant directly: every u<v produced by
  # .q_bd2_birth_sample() lies strictly inside a segment's interior --
  # confirmed both by construction (see .q_bd2_birth_sample's use of
  # segment-local composition sampling) and empirically below.
  set.seed(103)
  for (i in 1:2000) {
    rp <- .q_bd2_birth_sample(r, N, K_max = NULL, min_gap = 1L)
    added <- setdiff(which(rp == 1), which(r == 1))
    if (length(added) > 0) {
      new_pos <- added + 1L
      expect_false(any(new_pos %in% boundaries))
    }
  }
})

test_that("case 12: Type 3 birth is exactly normalized and matches 1/2 * 1/|E_B| * 1/|B_j|", {
  # Deterministic: enumerate EVERY target reachable from a small source
  # state by a single Type-3 birth (no enumeration inside the package
  # itself -- .q_bd2_birth_log_pmf() is closed-form -- this is just the
  # test's own brute-force reference), verify each against the exact
  # closed-form formula, and verify the birth-half probability mass
  # (including its own self-loop) sums to exactly 0.5.
  N <- 16
  r <- integer(N - 2)
  r[c(3, 12)] <- 1L
  cps <- 1 + which(r == 1)
  boundaries <- c(1L, cps, N)

  total_birth_mass <- 0
  n_targets <- 0
  for (seg in seq_len(length(boundaries) - 1L)) {
    lo <- boundaries[seg] + 1L; hi <- boundaries[seg + 1L] - 1L
    if (hi < lo + 1L) next
    for (u in lo:hi) for (v in lo:hi) {
      if (v <= u) next
      rt <- r; rt[c(u, v) - 1L] <- 1L
      if (!.r_is_admissible(rt, N)) next
      p <- exp(log_q_bd2_pmf(rt, r, N))
      total_birth_mass <- total_birth_mass + p
      n_targets <- n_targets + 1L
    }
  }
  expect_gt(n_targets, 0)
  expect_equal(total_birth_mass + 0.5 * .q_bd2_birth_self_mass(r, N, NULL, 1L), 0.5,
               tolerance = 1e-9)

  # light Monte Carlo as a loose secondary sanity check only (not the
  # primary correctness evidence -- that's the deterministic check above)
  set.seed(104)
  nsim <- 20000
  outs <- character(nsim)
  for (i in seq_len(nsim)) outs[i] <- paste(which(.q_bd2_core(r, N)$r_prop == 1), collapse = ",")
  tab <- table(outs) / nsim
  predicted <- sapply(names(tab), function(o) {
    idx <- if (nchar(o) == 0) integer(0) else as.integer(strsplit(o, ",")[[1]])
    rt <- integer(N - 2); rt[idx] <- 1L
    exp(log_q_bd2_pmf(rt, r, N))
  })
  expect_lt(max(abs(as.numeric(tab) - as.numeric(predicted))), 0.02)
})

test_that("case 13: Type 3 deletion is exactly normalized and matches 1/2 * 1/|D|", {
  # D(r) = all c(r)-1 adjacent-in-list pairs, unconditionally (verified:
  # deleting changepoints from an admissible state is always admissible --
  # see revision notes), so no admissibility filtering is needed to
  # construct the reference set here either.
  adjacent_pairs <- function(r) {
    ones <- which(r == 1)
    c_r <- length(ones)
    if (c_r < 2) return(list())
    lapply(seq_len(c_r - 1L), function(j) c(ones[j], ones[j + 1L]))
  }

  N <- 20
  r <- integer(N - 2)
  r[c(3, 6, 10, 14)] <- 1L # c(r) = 4, |D| = 3 adjacent pairs
  D <- adjacent_pairs(r)
  expect_equal(length(D), 3) # c(r) - 1

  total_death_mass <- 0
  for (pr in D) {
    rt <- r; rt[pr] <- 0L
    p <- exp(log_q_bd2_pmf(rt, r, N))
    expect_equal(p, 0.5 / length(D), tolerance = 1e-12) # exact, deterministic
    total_death_mass <- total_death_mass + p
  }
  p_death_self <- as.numeric(length(D) == 0)
  expect_equal(total_death_mass + 0.5 * p_death_self, 0.5, tolerance = 1e-9)
})

test_that("Type 3 (q_bd2) is exactly normalized including the diagonal, across several states", {
  # sum_{r'} q_bd2(r' | r) == 1 for every r' (all off-diagonal births,
  # all off-diagonal deaths, and the combined self-loop), deterministically.
  adjacent_pairs <- function(r) {
    ones <- which(r == 1)
    c_r <- length(ones)
    if (c_r < 2) return(list())
    lapply(seq_len(c_r - 1L), function(j) c(ones[j], ones[j + 1L]))
  }
  check_normalization <- function(r, N, K_max = NULL, min_gap = 1L) {
    cps <- 1 + which(r == 1)
    boundaries <- c(1L, cps, N)
    total <- 0
    for (seg in seq_len(length(boundaries) - 1L)) {
      lo <- boundaries[seg] + 1L; hi <- boundaries[seg + 1L] - 1L
      if (hi < lo + 1L) next
      for (u in lo:hi) for (v in lo:hi) {
        if (v <= u) next
        rt <- r; rt[c(u, v) - 1L] <- 1L
        if (!.r_is_admissible(rt, N, K_max, min_gap)) next
        total <- total + exp(log_q_bd2_pmf(rt, r, N, K_max, min_gap))
      }
    }
    D <- adjacent_pairs(r)
    for (pr in D) {
      rt <- r; rt[pr] <- 0L
      total <- total + exp(log_q_bd2_pmf(rt, r, N, K_max, min_gap))
    }
    total + exp(log_q_bd2_pmf(r, r, N, K_max, min_gap)) # + diagonal
  }

  N <- 18
  expect_equal(check_normalization(integer(N - 2), N), 1, tolerance = 1e-9)
  r1 <- integer(N - 2); r1[c(4, 9, 13)] <- 1L
  expect_equal(check_normalization(r1, N), 1, tolerance = 1e-9)
  r2 <- integer(N - 2); r2[c(3, 4, 5)] <- 1L # tightly clustered
  expect_equal(check_normalization(r2, N), 1, tolerance = 1e-9)
  r3 <- integer(N - 2); r3[1:10] <- 1L # near-saturated
  expect_equal(check_normalization(r3, N), 1, tolerance = 1e-9)
  expect_equal(check_normalization(r1, N, K_max = 4, min_gap = 2L), 1, tolerance = 1e-9)
})

test_that("case 14: singleton candidate sets are sampled correctly (no sample() gotcha)", {
  # sample(x, 1) on a length-1 numeric x treats x as "sample from 1:x", not
  # "return the literal value x" -- a classic base-R gotcha. Every propose
  # function in this package uses sample.int(length(x), 1) + indexing
  # instead. This test constructs states where a candidate set has exactly
  # one element and confirms the singleton is always returned deterministically.
  set.seed(106)
  N <- 20

  # Type 2 death: c(r) = 1 -> exactly one changepoint to delete. The
  # log-pmf's self-loop/death mass should be exactly 0.5 (death always
  # succeeds deterministically on the single candidate), not distorted by
  # sample()'s length-1 gotcha.
  r <- integer(N - 2); r[7] <- 1L
  expect_equal(exp(log_q_bd_pmf(integer(N - 2), r, N)), 0.5, tolerance = 1e-8)

  # Type 3 death: c(r) = 2 -> exactly one adjacent pair (|D| = c-1 = 1)
  r3 <- integer(N - 2); r3[c(5, 6)] <- 1L
  expect_equal(sum(r3) - 1L, 1) # |D| = c-1
  target_after_death <- integer(N - 2)
  expect_equal(exp(log_q_bd2_pmf(target_after_death, r3, N)), 0.5, tolerance = 1e-8)

  # confirm empirically over many draws that the single candidate is what's
  # actually produced whenever the death branch is taken
  n_death_events <- 0L
  for (i in 1:5000) {
    rp <- .q_bd2_core(r3, N)$r_prop
    if (sum(rp) == 0) n_death_events <- n_death_events + 1L
    expect_true(sum(rp) %in% c(0, 2, 4)) # death(0), self-loop(2), or birth(4)
  }
  expect_gt(n_death_events, 0)
})

test_that("case 15: Type 4 pointwise probability agrees with 1/[c(r)(n-2-c(r))]", {
  # Deterministic primary check (per Item F): enumerate every possible
  # single-swap target from a small source state directly, rather than
  # relying on Monte Carlo sampling frequencies for the core correctness
  # claim.
  N <- 14
  r <- integer(N - 2)
  r[c(3, 8)] <- 1L
  c_r <- sum(r)
  ones <- which(r == 1)
  zeros <- which(r == 0)
  expected_p <- 1 / (c_r * (N - 2 - c_r))

  total_off_diag <- 0
  n_targets <- 0
  for (s in ones) {
    for (ss in zeros) {
      rt <- r
      rt[s] <- 0L
      rt[ss] <- 1L
      if (identical(rt, r)) next # shouldn't happen (s != ss by construction)
      p <- exp(log_q_shift_pmf(rt, r, N))
      if (.r_is_admissible(rt, N)) {
        expect_equal(p, expected_p, tolerance = 1e-12) # exact, deterministic
        total_off_diag <- total_off_diag + p
        n_targets <- n_targets + 1L
      } else {
        expect_equal(p, 0, tolerance = 1e-12)
      }
    }
  }
  expect_gt(n_targets, 0)
  self_p <- exp(log_q_shift_pmf(r, r, N))
  expect_equal(total_off_diag + self_p, 1, tolerance = 1e-9) # full normalization

  # light Monte Carlo as a loose secondary sanity check only
  set.seed(107)
  nsim <- 20000
  outs <- character(nsim)
  for (i in seq_len(nsim)) outs[i] <- paste(which(q_shift(r, N) == 1), collapse = ",")
  tab <- table(outs) / nsim
  non_self <- tab[names(tab) != paste(which(r == 1), collapse = ",")]
  expect_true(all(abs(as.numeric(non_self) - expected_p) < 0.01)) # loose, secondary
})

# ============================================================================
# Detailed balance: exhaustive verification on a small enumerated state
# space, for each individual kernel and for the fixed mixture, under an
# arbitrary positive toy target (cases 16-17).
# ============================================================================

.enumerate_admissible_states <- function(N, K_max = NULL, min_gap = 1L) {
  M <- N - 2L
  states <- list()
  for (code in 0:(2^M - 1)) {
    r <- as.integer(intToBits(code))[1:M]
    if (.r_is_admissible(r, N, K_max, min_gap)) states[[length(states) + 1]] <- r
  }
  states
}

.toy_target <- function(r) {
  exp(0.7 * sum(r) - 0.05 * sum((which(r == 1))^1.3) + 0.3 * sin(sum(r) * 1.7 + 1))
}

test_that("case 16: each individual proposal kernel satisfies detailed balance exactly", {
  N <- 8 # small enough to enumerate exhaustively (M = 6, up to 64 raw states)
  states <- .enumerate_admissible_states(N)
  pi_vals <- sapply(states, .toy_target)
  ns <- length(states)
  expect_gt(ns, 5) # sanity: nontrivial state space

  for (type in 1:4) {
    max_violation <- 0
    for (i in seq_len(ns)) {
      for (j in seq_len(ns)) {
        if (i == j) next
        qij <- exp(pproposal(type, states[[j]], states[[i]], N, 1 / 30, 0.05))
        qji <- exp(pproposal(type, states[[i]], states[[j]], N, 1 / 30, 0.05))
        if (qij == 0 && qji == 0) next
        alpha_ij <- if (qij == 0) 0 else min(1, (pi_vals[j] * qji) / (pi_vals[i] * qij))
        alpha_ji <- if (qji == 0) 0 else min(1, (pi_vals[i] * qij) / (pi_vals[j] * qji))
        lhs <- pi_vals[i] * qij * alpha_ij
        rhs <- pi_vals[j] * qji * alpha_ji
        max_violation <- max(max_violation, abs(lhs - rhs))
      }
    }
    expect_lt(max_violation, 1e-9)
  }
})

test_that("case 17: the fixed-weight mixture kernel satisfies detailed balance exactly", {
  N <- 8
  states <- .enumerate_admissible_states(N)
  pi_vals <- sapply(states, .toy_target)
  ns <- length(states)
  weights <- c(1 / 4, 1 / 8, 1 / 8, 1 / 2)

  max_violation <- 0
  for (i in seq_len(ns)) {
    for (j in seq_len(ns)) {
      if (i == j) next
      Pij <- 0; Pji <- 0
      for (type in 1:4) {
        qij <- exp(pproposal(type, states[[j]], states[[i]], N, 1 / 30, 0.05))
        qji <- exp(pproposal(type, states[[i]], states[[j]], N, 1 / 30, 0.05))
        alpha_ij <- if (qij == 0) 0 else min(1, (pi_vals[j] * qji) / (pi_vals[i] * qij))
        alpha_ji <- if (qji == 0) 0 else min(1, (pi_vals[i] * qij) / (pi_vals[j] * qji))
        Pij <- Pij + weights[type] * qij * alpha_ij
        Pji <- Pji + weights[type] * qji * alpha_ji
      }
      max_violation <- max(max_violation, abs(pi_vals[i] * Pij - pi_vals[j] * Pji))
    }
  }
  expect_lt(max_violation, 1e-9)
})

test_that("detailed balance also holds with K_max and min_gap binding", {
  N <- 9
  states <- .enumerate_admissible_states(N, K_max = 4, min_gap = 2L)
  pi_vals <- sapply(states, .toy_target)
  ns <- length(states)
  expect_gt(ns, 3)
  weights <- c(1 / 4, 1 / 8, 1 / 8, 1 / 2)

  max_violation <- 0
  for (i in seq_len(ns)) {
    for (j in seq_len(ns)) {
      if (i == j) next
      Pij <- 0; Pji <- 0
      for (type in 1:4) {
        qij <- exp(pproposal(type, states[[j]], states[[i]], N, 1 / 30, 0.05, K_max = 4, min_gap = 2L))
        qji <- exp(pproposal(type, states[[i]], states[[j]], N, 1 / 30, 0.05, K_max = 4, min_gap = 2L))
        alpha_ij <- if (qij == 0) 0 else min(1, (pi_vals[j] * qji) / (pi_vals[i] * qij))
        alpha_ji <- if (qji == 0) 0 else min(1, (pi_vals[i] * qij) / (pi_vals[j] * qji))
        Pij <- Pij + weights[type] * qij * alpha_ij
        Pji <- Pji + weights[type] * qji * alpha_ji
      }
      max_violation <- max(max_violation, abs(pi_vals[i] * Pij - pi_vals[j] * Pji))
    }
  }
  expect_lt(max_violation, 1e-9)
})

test_that("MHsearch runs end-to-end with the new proposal mechanism", {
  set.seed(108)
  n <- 80
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  x <- cumsum(rnorm(n, 0, 0.02))
  y <- cumsum(rnorm(n, 0, 0.02))

  res <- MHsearch(t, x, y, dt = 0.05, iter_max = 300, burn_in = 50, show_progress = FALSE)
  expect_true(length(res$cps_list) > 0)
  expect_true(all(is.finite(res$info_table$score_cur)))

  # with an explicit K_max/min_gap supplied end to end
  res2 <- MHsearch(t, x, y, dt = 0.05, iter_max = 300, burn_in = 50, show_progress = FALSE,
                    K_max = 6, min_gap = 2L)
  expect_true(length(res2$cps_list) > 0)
  n_segs <- vapply(res2$cps_list, function(z) length(z$this_cp) + 1L, integer(1))
  expect_true(all(n_segs <= 6))
})

# ============================================================================
# Backward-compatibility wrappers (Item A)
# ============================================================================

test_that("q_new(lambda_r, dt, N) restores the original 3-arg calling convention", {
  set.seed(201)
  r <- q_new(1 / 30, 0.05, 20)
  expect_length(r, 18)
  expect_true(all(r %in% c(0L, 1L)))
})

test_that("q_bd(r_cur) and q_bd2(r_cur) restore the original list(r_prop=, status=) structure", {
  set.seed(202)
  N <- 20
  r_cur <- integer(N - 2); r_cur[c(5, 10)] <- 1L

  res <- q_bd(r_cur)
  expect_named(res, c("r_prop", "status"))
  expect_length(res$r_prop, N - 2)
  expect_true(res$status %in% c(0, 1))

  res2 <- q_bd2(r_cur)
  expect_named(res2, c("r_prop", "status"))
  expect_length(res2$r_prop, N - 2)
  expect_true(res2$status %in% c(0, 1))
})

test_that("q_shift(r_cur) keeps returning r_prop directly (unchanged structure)", {
  set.seed(203)
  N <- 20
  r_cur <- integer(N - 2); r_cur[c(5, 10)] <- 1L
  r <- q_shift(r_cur)
  expect_length(r, N - 2)
  expect_true(is.numeric(r))
})

test_that("q_as(r_cur)/q_ds(r_cur) are restored, dispatching to corrected Type-3 machinery", {
  set.seed(204)
  N <- 20
  r_cur <- integer(N - 2); r_cur[c(5, 10)] <- 1L

  r_birth <- q_as(r_cur)
  expect_length(r_birth, N - 2)
  # q_as always attempts a birth: result has >= as many changepoints
  expect_gte(sum(r_birth), sum(r_cur))

  r_death <- q_ds(r_cur)
  expect_length(r_death, N - 2)
  expect_lte(sum(r_death), sum(r_cur))
})

test_that("q_as_pmf() correctly requires (r_target, r_source), unlike the original 1-arg form", {
  N <- 20
  r_cur <- integer(N - 2); r_cur[c(5, 10)] <- 1L
  target <- r_cur; target[c(7, 8)] <- 1L
  p <- q_as_pmf(target, r_cur)
  expect_true(is.numeric(p) && p >= 0 && p <= 1)
  # matches the corrected internal birth density exactly (not a log)
  expect_equal(p, exp(.q_bd2_birth_log_pmf(target, r_cur, N, NULL, 1L)))
})

test_that("q_ds_pmf(r) is fully restored with the original 1-argument signature, now correct", {
  N <- 20
  r <- integer(N - 2); r[c(5, 6, 10)] <- 1L # |D| = 2
  expect_equal(q_ds_pmf(r), 0.5 / 2)
  expect_equal(q_ds_pmf(integer(N - 2)), 0) # no changepoints -> no legal deletion
  r_one <- integer(N - 2); r_one[5] <- 1L
  expect_equal(q_ds_pmf(r_one), 0) # only 1 changepoint -> no legal deletion
})

test_that("compat wrappers reach the SAME corrected machinery detailed balance already covers", {
  # q_as()/q_ds() are literally the birth-only/death-only halves of
  # .q_bd2_core(), which detailed-balance tests (16-17, above) already
  # cover as part of the full Type-3 mixture -- this just confirms the
  # wrapper dispatch itself is wired correctly, via a direct equivalence
  # check against the internal functions they're documented to call.
  N <- 16
  r <- integer(N - 2); r[c(3, 12)] <- 1L
  set.seed(205)
  s1 <- .Random.seed
  r1 <- q_as(r, N)
  s2 <- .Random.seed
  set.seed(205)
  r2 <- .q_bd2_birth_sample(r, N, NULL, 1L)
  expect_identical(r1, r2)
})

test_that("q_new(lambda_r, dt, N) retains the TRUE original raw-generator behavior", {
  # Length and marginal distribution should match an unconditional
  # Bernoulli(N-2, p) draw exactly -- no admissibility filtering, no
  # self-transition, no dependence on any current state.
  set.seed(301)
  N <- 20
  lambda_r <- 1 / 30
  dt <- 0.05
  p <- 1 - exp(-lambda_r * dt)

  r <- q_new(lambda_r, dt, N)
  expect_length(r, N - 2)
  expect_true(all(r %in% c(0L, 1L)))

  # empirical per-position marginal probability matches p (loose sanity
  # check; this is a plain iid Bernoulli generator, nothing more)
  nsim <- 20000
  draws <- replicate(nsim, q_new(lambda_r, dt, N))
  expect_equal(mean(draws), p, tolerance = 0.01)

  # crucially: it can (and, at these settings, essentially always will)
  # produce the all-zero vector when the raw draw happens to be all-zero,
  # WITHOUT reference to any "current state" -- there is no such concept
  # in this wrapper's contract. Confirm it does NOT go through the
  # admissibility-aware .q_new_core() by checking it can return a raw
  # draw that violates c(r) <= n-3 (the saturated / all-ones case),
  # which the corrected kernel would never emit as-is.
  set.seed(302)
  saturated_seen <- FALSE
  for (i in 1:200) {
    r_hi <- q_new(lambda_r = 50, dt = 1, N) # p essentially 1
    if (sum(r_hi) == N - 2) { saturated_seen <- TRUE; break }
  }
  expect_true(saturated_seen) # raw generator: no admissibility filtering
})

test_that("Type 3 death simplification: exact for c-1 adjacent pairs, with K_max/min_gap binding", {
  # Item 2: D(r) = all c-1 adjacent pairs unconditionally; verify the
  # closed-form death density and normalization hold exactly even when
  # K_max/min_gap are binding elsewhere in the state (they cannot affect
  # the death branch's admissibility per the proof, but confirm the
  # THIS-BRANCH VALUE is correct end-to-end anyway).
  N <- 22
  r <- integer(N - 2)
  r[c(3, 6, 9, 13, 17)] <- 1L # c(r) = 5, |D| = 4
  c_r <- sum(r)

  for (K_max in list(NULL, 6, 8)) {
    for (min_gap in c(1L, 2L, 3L)) {
      ones <- which(r == 1)
      total <- 0
      for (j in seq_len(c_r - 1L)) {
        rt <- r; rt[c(ones[j], ones[j + 1L])] <- 0L
        p <- exp(log_q_bd2_pmf(rt, r, N, K_max, min_gap))
        expect_equal(p, 0.5 / (c_r - 1L), tolerance = 1e-12)
        total <- total + p
      }
      expect_equal(total, 0.5, tolerance = 1e-9) # death half fully accounted for
    }
  }

  # c < 2: death half is a pure self-transition (mass 0.5)
  r_one <- integer(N - 2); r_one[5] <- 1L
  expect_equal(exp(log_q_bd2_pmf(r_one, r_one, N)) -
                 0.5 * .q_bd2_birth_self_mass(r_one, N, NULL, 1L), 0.5, tolerance = 1e-9)

  # a non-adjacent-in-list pair must NOT be treated as a legal Type-3 death
  r2 <- integer(N - 2); r2[c(3, 6, 9)] <- 1L # ones at 3,6,9
  bad_target <- r2; bad_target[c(3, 9)] <- 0L # removes {3,9}, NOT adjacent (6 is between)
  expect_equal(exp(log_q_bd2_pmf(bad_target, r2, N)), 0)
})
