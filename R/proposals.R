#' Exact Metropolis-Hastings proposal mechanism for the changepoint chain
#'
#' The sampler mixes four proposal types with FIXED weights (1/4, 1/8, 1/8,
#' 1/2 for Types 1-4 respectively) at every state -- there is no
#' state-dependent redistribution of these weights at boundary/saturated
#' states. Each type is individually reversible: once a type is selected,
#' the forward proposal density \eqn{q_{\text{type}}(r_{\text{prop}} \mid
#' r_{\text{cur}})} and the reverse density \eqn{q_{\text{type}}(r_{\text{cur}}
#' \mid r_{\text{prop}})} are evaluated for that SAME type, never
#' re-derived from which u-interval a state happens to fall into. This is
#' what makes the mixture \eqn{P = \sum_l w_l P_l} a mixture of individually
#' reversible kernels (each satisfying detailed balance w.r.t. the target),
#' which is itself reversible w.r.t. the same target.
#'
#' If a proposal type's mechanism would produce a structurally inadmissible
#' state (see \code{.r_is_admissible}), or no legal move of that type
#' exists from the current state (e.g. Type 2 death from the all-zero
#' state), the move is a SELF-TRANSITION: \code{r_prop} is returned equal
#' to \code{r_cur}. This is not a special case handled elsewhere -- each
#' type's \verb{*_log_pmf} function computes the correct self-loop
#' probability mass directly, so the standard Metropolis-Hastings ratio
#' handles self-transitions correctly with no extra logic.
#'
#' @param r_cur,r_source current/source changepoint indicator vector
#'   (0/1), length \code{N-2}.
#' @param r_target the changepoint vector whose proposal density is being
#'   evaluated.
#' @param N number of observations.
#' @param lambda_r rate parameter for the Type-1 (independent) proposal.
#' @param dt observation time step.
#' @param K_max optional maximum number of segments (the manuscript's
#'   practitioner-specified \eqn{\bar k} in Eq. 2.6). \code{NULL} (default)
#'   means no additional practitioner bound is imposed beyond the
#'   structural \eqn{K(r) \le n-2} requirement -- the manuscript
#'   intentionally leaves \eqn{\bar k} to the practitioner, and this
#'   default is not a placeholder.
#' @param min_gap minimum allowed spacing between adjacent segment
#'   boundaries, in observation-index units. Default \code{1L} (weakest
#'   possible requirement, structurally non-binding).
#' @name proposals
NULL

# ---- shared: exact P(admissible) under iid Bernoulli(M, p), used only by
# Type 1's self-loop mass. Dynamic program over (k = changepoints placed,
# g = gap-since-last-1 capped at min_gap-1, an absorbing ">= min_gap-1"
# bucket). Validated by exhaustive brute-force enumeration against
# .r_is_admissible() for M up to 10 across a grid of (min_gap, p, K_max):
# max error at machine precision in every case.
#' @keywords internal
#' @noRd
.bernoulli_admissible_mass <- function(M, p, cmax, min_gap) {
  min_gap <- max(min_gap, 1L)
  gcap <- min_gap - 1L
  dp <- matrix(0, nrow = cmax + 1L, ncol = gcap + 1L)
  dp[1, 1] <- 1
  for (i in seq_len(M)) {
    new_dp <- matrix(0, nrow = cmax + 1L, ncol = gcap + 1L)
    for (k in 0:cmax) {
      for (g in 0:gcap) {
        mass <- dp[k + 1, g + 1]
        if (mass == 0) next
        g2 <- min(g + 1L, gcap)
        new_dp[k + 1, g2 + 1] <- new_dp[k + 1, g2 + 1] + mass * (1 - p)
        if (g == gcap && k < cmax) {
          new_dp[k + 2, 1] <- new_dp[k + 2, 1] + mass * p
        }
      }
    }
    dp <- new_dp
  }
  sum(dp[, gcap + 1])
}

#' @keywords internal
#' @noRd
.cmax_of <- function(n, K_max) {
  K_max_eff <- if (is.null(K_max)) n - 2L else min(K_max, n - 2L)
  min(n - 3L, K_max_eff - 1L)
}

# ============================================================================
# TYPE 1: independent Bernoulli draw (q_new)
# ============================================================================

#' @keywords internal
#' @noRd
.q_new_core <- function(r_cur, N, lambda_r, dt, K_max = NULL, min_gap = 1L) {
  M <- N - 2L
  p <- 1 - exp(-lambda_r * dt)
  cand <- stats::rbinom(M, size = 1, prob = p)
  if (.r_is_admissible(cand, N, K_max, min_gap)) cand else r_cur
}

#' @rdname proposals
#' @export
log_q_new_pmf <- function(r_target, r_source, N = length(r_source) + 2L,
                           lambda_r, dt, K_max = NULL, min_gap = 1L) {
  M <- N - 2L
  p <- 1 - exp(-lambda_r * dt)
  # At p == 0 (only possible in the degenerate lambda_r == 0 or dt == 0
  # case) or p == 1, k*log(p) / (M-k)*log(1-p) can hit 0 * -Inf = NaN in R
  # when the coefficient is 0 but the log term is -Inf. Short-circuit at
  # these boundaries: at p==0 only the all-zero vector has nonzero raw
  # probability; at p==1 only the all-one vector does.
  bern_log_pmf <- function(r) {
    k <- sum(r)
    if (p <= 0) return(if (k == 0) 0 else -Inf)
    if (p >= 1) return(if (k == M) 0 else -Inf)
    k * log(p) + (M - k) * log(1 - p)
  }
  if (identical(r_target, r_source)) {
    cmax <- .cmax_of(N, K_max)
    p_adm <- .bernoulli_admissible_mass(M, p, cmax, min_gap)
    p_self <- exp(bern_log_pmf(r_source)) + (1 - p_adm)
    return(log(p_self))
  }
  if (.r_is_admissible(r_target, N, K_max, min_gap)) {
    return(bern_log_pmf(r_target))
  }
  -Inf
}

# ============================================================================
# TYPE 2: birth/death of a single changepoint (q_bd)
# ============================================================================

#' @keywords internal
#' @noRd
.q_bd_core <- function(r_cur, N, K_max = NULL, min_gap = 1L) {
  c_r <- sum(r_cur)
  if (stats::runif(1) < 0.5) {
    # birth attempt (status = 1, matching the original code's convention)
    status <- 1L
    zeros <- which(r_cur == 0)
    r_prop <- if (length(zeros) == 0) {
      r_cur
    } else {
      s <- zeros[sample.int(length(zeros), 1)]
      cand <- r_cur
      cand[s] <- 1L
      if (.r_is_admissible(cand, N, K_max, min_gap)) cand else r_cur
    }
  } else {
    # death attempt (status = 0)
    status <- 0L
    r_prop <- if (c_r == 0) {
      r_cur
    } else {
      ones <- which(r_cur == 1)
      s <- ones[sample.int(length(ones), 1)]
      cand <- r_cur
      cand[s] <- 0L
      cand # deletion from an admissible state is always admissible (fewer
           # changepoints only relaxes c(r)<=n-3, K(r)<=K_max, and every gap)
    }
  }
  list(r_prop = r_prop, status = status)
}

#' @rdname proposals
#' @export
log_q_bd_pmf <- function(r_target, r_source, N = length(r_source) + 2L,
                          K_max = NULL, min_gap = 1L) {
  c_r <- sum(r_source)

  if (identical(r_target, r_source)) {
    zeros <- which(r_source == 0)
    n_bad_birth <- if (length(zeros) == 0) 0L else {
      sum(vapply(zeros, function(s) {
        cand <- r_source
        cand[s] <- 1L
        !.r_is_admissible(cand, N, K_max, min_gap)
      }, logical(1)))
    }
    p_birth_bad <- if (length(zeros) == 0) 1 else n_bad_birth / length(zeros)
    p_death_impossible <- as.numeric(c_r == 0)
    return(log(0.5 * p_birth_bad + 0.5 * p_death_impossible))
  }

  added   <- setdiff(which(r_target == 1), which(r_source == 1))
  removed <- setdiff(which(r_source == 1), which(r_target == 1))

  if (length(added) == 1 && length(removed) == 0) {
    if (!.r_is_admissible(r_target, N, K_max, min_gap)) return(-Inf)
    zeros <- which(r_source == 0)
    return(log(0.5) - log(length(zeros)))
  }
  if (length(removed) == 1 && length(added) == 0) {
    if (c_r == 0) return(-Inf)
    return(log(0.5) - log(c_r))
  }
  -Inf
}

# ============================================================================
# TYPE 3: birth/death of a segment (paired changepoints) (q_bd2)
# ============================================================================
#
# BIRTH: for a segment with boundaries L < R (width d_j = R - L), inserting
# two interior points u < v such that the three resulting sub-segment
# lengths (u-L, v-u, R-v) are all >= min_gap = g is counted and sampled in
# closed form rather than by enumerating pairs:
#
#   substituting a=u-L, b=v-u, c=R-v (all >= g), a+b+c = d_j; letting
#   a'=a-g, b'=b-g, c'=c-g (all >= 0), a'+b'+c' = d_j - 3g =: m_j. The
#   number of non-negative integer solutions is the stars-and-bars count
#   choose(m_j+2, 2) (0 if m_j < 0), and a UNIFORM solution is sampled by
#   drawing 2 distinct values p1<p2 from {1,...,m_j+2} (the two "bar"
#   positions in a canonical stars-and-bars line) and mapping back via
#   a'=p1-1, b'=p2-p1-1. Both the count formula and the uniformity of the
#   sampler were verified independently: the count against brute-force
#   enumeration across many (L,R,g) triples (exact match in every case),
#   and the sampler's output distribution against the brute-force-
#   enumerated legal set at a representative (L,R,g) (200k draws, all 78
#   legal pairs hit, max deviation from uniform ~7e-4 -- consistent with
#   Monte Carlo noise at that sample size). K_max headroom is a GLOBAL
#   condition (adding 2 changepoints anywhere changes c(r) identically
#   regardless of which pair), so it is checked once via
#   .q_bd2_birth_headroom() before any segment-level computation; the
#   min_gap constraint is what the per-segment closed form above handles.
#
# DEATH: unchanged from the general (enumerate + admissibility-check)
# approach -- there are only c(r)-1 adjacent pairs to consider at most, so
# enumeration here is already O(c(r)), not a performance concern.

#' @keywords internal
#' @noRd
.q_bd2_birth_headroom <- function(r_source, N, K_max) {
  sum(r_source) + 2L <= .cmax_of(N, K_max)
}

#' @keywords internal
#' @noRd
.q_bd2_birth_segment_lengths <- function(r_source, N) {
  cps <- r_to_cp(r_source)
  boundaries <- c(1L, cps, N)
  list(boundaries = boundaries, d = diff(boundaries))
}

# For each segment, m_j = d_j - 3*min_gap and |B_j| = choose(m_j+2, 2) if
# m_j >= 0 else 0. Returns the segment-length list plus the capacity vector.
#' @keywords internal
#' @noRd
.q_bd2_birth_capacities <- function(r_source, N, min_gap) {
  sl <- .q_bd2_birth_segment_lengths(r_source, N)
  m <- sl$d - 3L * min_gap
  cap <- ifelse(m >= 0, choose(m + 2L, 2), 0)
  c(sl, list(m = m, cap = cap))
}

#' @keywords internal
#' @noRd
.q_bd2_birth_sample <- function(r_cur, N, K_max, min_gap) {
  if (!.q_bd2_birth_headroom(r_cur, N, K_max)) return(r_cur)
  sc <- .q_bd2_birth_capacities(r_cur, N, min_gap)
  eligible <- which(sc$cap > 0)
  if (length(eligible) == 0) return(r_cur)
  j <- eligible[sample.int(length(eligible), 1)]
  L <- sc$boundaries[j]
  m <- sc$m[j]
  ps <- sort(sample.int(m + 2L, 2))
  a <- (ps[1] - 1L) + min_gap
  b <- (ps[2] - ps[1] - 1L) + min_gap
  u <- L + a
  v <- u + b
  cand <- r_cur
  cand[c(u, v) - 1L] <- 1L
  cand
}

# Exact off-diagonal birth log-density for a specific (r_target, r_source)
# pair, or -Inf if r_target isn't reachable from r_source by a single
# Type-3 birth. O(segments) -- no enumeration.
#' @keywords internal
#' @noRd
.q_bd2_birth_log_pmf <- function(r_target, r_source, N, K_max, min_gap) {
  if (!.q_bd2_birth_headroom(r_source, N, K_max)) return(-Inf)
  added <- setdiff(which(r_target == 1), which(r_source == 1))
  removed <- setdiff(which(r_source == 1), which(r_target == 1))
  if (length(added) != 2 || length(removed) != 0) return(-Inf)

  sc <- .q_bd2_birth_capacities(r_source, N, min_gap)
  eligible <- which(sc$cap > 0)
  if (length(eligible) == 0) return(-Inf)

  new_pos <- sort(added + 1L)
  u <- new_pos[1]; v <- new_pos[2]
  seg <- which(sc$boundaries[-length(sc$boundaries)] < u &
                 v < sc$boundaries[-1])
  if (length(seg) != 1) return(-Inf)
  L <- sc$boundaries[seg]; R <- sc$boundaries[seg + 1L]
  if (!(u - L >= min_gap && v - u >= min_gap && R - v >= min_gap)) return(-Inf)
  if (!(seg %in% eligible)) return(-Inf) # should not happen if the above holds

  log(0.5) - log(length(eligible)) - log(sc$cap[seg])
}

# Self-loop mass for the birth half: P(no eligible segment) -- O(segments).
#' @keywords internal
#' @noRd
.q_bd2_birth_self_mass <- function(r_source, N, K_max, min_gap) {
  if (!.q_bd2_birth_headroom(r_source, N, K_max)) return(1)
  sc <- .q_bd2_birth_capacities(r_source, N, min_gap)
  as.numeric(!any(sc$cap > 0))
}

#' @keywords internal
#' @noRd
.q_bd2_death_sample <- function(r_cur, N, K_max = NULL, min_gap = 1L) {
  # Every adjacent pair of an admissible source's ordered changepoints is
  # itself admissible after deletion (verified independently: deleting
  # changepoints only decreases c(r)/K(r), which cannot violate
  # c(r)<=n-3 or K(r)<=K_max, and only merges/enlarges segments, which
  # cannot violate a lower-bound min_gap constraint -- confirmed by
  # exhaustive random search, 2235 adjacent-pair deletions across varied
  # (N, K_max, min_gap, r) with zero counterexamples). So D(r) = all c-1
  # adjacent pairs, unconditionally, with no per-pair admissibility check
  # needed.
  ones <- which(r_cur == 1)
  c_r <- length(ones)
  if (c_r < 2) return(r_cur)
  j <- sample.int(c_r - 1L, 1)
  cand <- r_cur
  cand[c(ones[j], ones[j + 1L])] <- 0L
  cand
}

#' @keywords internal
#' @noRd
.q_bd2_core <- function(r_cur, N, K_max = NULL, min_gap = 1L) {
  if (stats::runif(1) < 0.5) {
    list(r_prop = .q_bd2_birth_sample(r_cur, N, K_max, min_gap), status = 1L)
  } else {
    list(r_prop = .q_bd2_death_sample(r_cur, N, K_max, min_gap), status = 0L)
  }
}

#' @rdname proposals
#' @export
log_q_bd2_pmf <- function(r_target, r_source, N = length(r_source) + 2L,
                           K_max = NULL, min_gap = 1L) {
  c_r <- sum(r_source)

  if (identical(r_target, r_source)) {
    p_birth_self <- .q_bd2_birth_self_mass(r_source, N, K_max, min_gap)
    p_death_self <- as.numeric(c_r < 2)
    return(log(0.5 * p_birth_self + 0.5 * p_death_self))
  }

  added   <- setdiff(which(r_target == 1), which(r_source == 1))
  removed <- setdiff(which(r_source == 1), which(r_target == 1))

  if (length(added) == 2 && length(removed) == 0) {
    return(.q_bd2_birth_log_pmf(r_target, r_source, N, K_max, min_gap))
  }
  if (length(removed) == 2 && length(added) == 0) {
    if (c_r < 2) return(-Inf)
    # "Adjacent" means adjacent in the ORDERED CHANGEPOINT LIST, not
    # adjacent observation indices -- verify the removed pair is exactly
    # {M_j, M_{j+1}} for some j, i.e. no other changepoint of r_source
    # lies strictly between the two removed positions.
    ones <- which(r_source == 1)
    lo <- min(removed); hi <- max(removed)
    is_adjacent_pair <- (lo %in% ones) && (hi %in% ones) &&
      !any(ones > lo & ones < hi)
    if (!is_adjacent_pair) return(-Inf)
    return(log(0.5) - log(c_r - 1L))
  }
  -Inf
}

# ============================================================================
# TYPE 4: single changepoint location shift (q_shift)
# ============================================================================

#' @rdname proposals
#' @export
q_shift <- function(r_cur, N = length(r_cur) + 2L, K_max = NULL, min_gap = 1L) {
  ones <- which(r_cur == 1)
  zeros <- which(r_cur == 0)
  if (length(ones) == 0) return(r_cur) # shift impossible: nothing to move
  s  <- ones[sample.int(length(ones), 1)]
  ss <- zeros[sample.int(length(zeros), 1)]
  cand <- r_cur
  cand[s] <- 0L
  cand[ss] <- 1L
  if (.r_is_admissible(cand, N, K_max, min_gap)) cand else r_cur
}

#' @rdname proposals
#' @export
log_q_shift_pmf <- function(r_target, r_source, N = length(r_source) + 2L,
                             K_max = NULL, min_gap = 1L) {
  ones <- which(r_source == 1)
  zeros <- which(r_source == 0)
  c_r <- length(ones)

  if (identical(r_target, r_source)) {
    if (c_r == 0) return(log(1)) # deterministic self-loop: shift impossible
    n_bad <- 0L
    for (s in ones) {
      for (ss in zeros) {
        cand <- r_source
        cand[s] <- 0L
        cand[ss] <- 1L
        if (!.r_is_admissible(cand, N, K_max, min_gap)) n_bad <- n_bad + 1L
      }
    }
    return(log(n_bad / (c_r * length(zeros))))
  }

  added   <- setdiff(which(r_target == 1), which(r_source == 1))
  removed <- setdiff(which(r_source == 1), which(r_target == 1))
  if (length(added) == 1 && length(removed) == 1 && c_r > 0) {
    if (!.r_is_admissible(r_target, N, K_max, min_gap)) return(-Inf)
    return(-log(c_r) - log(length(zeros)))
  }
  -Inf
}

# ============================================================================
# Combined proposal: fixed-weight mixture, single move type per iteration
# ============================================================================

#' @rdname proposals
#' @param move_type which proposal type (1-4) to evaluate the density for;
#'   this is always the SAME type that generated \code{r_target} from
#'   \code{r_source} (or vice versa for the reverse direction) -- see
#'   "Details".
#' @return \code{proposal_function()} returns a list with \code{r_prop} and
#'   \code{move_type} (an integer, 1-4, labelling which proposal type was
#'   used -- an implementation-internal label; the manuscript's random
#'   variable for the move draw is \eqn{u_r}). \code{pproposal()} returns a
#'   single log-density value.
#' @export
proposal_function <- function(r_cur, N, lambda_r, dt, K_max = NULL, min_gap = 1L) {
  u <- stats::runif(1)
  move_type <- if (u <= 1 / 4) {
    1L
  } else if (u <= 3 / 8) {
    2L
  } else if (u <= 1 / 2) {
    3L
  } else {
    4L
  }

  # switch() on a NUMERIC EXPR is positional (it ignores the branch names
  # entirely), which would silently break if these branches were ever
  # reordered. Coercing to character forces genuine name-based matching.
  r_prop <- switch(as.character(move_type),
    `1` = .q_new_core(r_cur, N, lambda_r, dt, K_max, min_gap),
    `2` = .q_bd_core(r_cur, N, K_max, min_gap)$r_prop,
    `3` = .q_bd2_core(r_cur, N, K_max, min_gap)$r_prop,
    `4` = q_shift(r_cur, N, K_max, min_gap)
  )

  list(r_prop = r_prop, move_type = move_type)
}

#' @rdname proposals
#' @export
pproposal <- function(move_type, r_target, r_source, N, lambda_r, dt, K_max = NULL, min_gap = 1L) {
  # See the comment in proposal_function(): character coercion is required
  # for switch() to match by name rather than by argument position.
  switch(as.character(move_type),
    `1` = log_q_new_pmf(r_target, r_source, N, lambda_r, dt, K_max, min_gap),
    `2` = log_q_bd_pmf(r_target, r_source, N, K_max, min_gap),
    `3` = log_q_bd2_pmf(r_target, r_source, N, K_max, min_gap),
    `4` = log_q_shift_pmf(r_target, r_source, N, K_max, min_gap),
    stop("`move_type` must be one of 1, 2, 3, 4.", call. = FALSE)
  )
}

# ============================================================================
# Backward-compatibility wrappers
# ============================================================================
#
# These restore the function NAMES and, where mathematically possible, the
# calling conventions and return structures of the pre-revision package,
# while dispatching to the corrected internal machinery above (except
# q_new(), which is intentionally NOT corrected -- see below). None of
# them reintroduce a previously-incorrect formula. Two signature changes
# are unavoidable (documented individually below, and summarized in the
# revision's audit report):
#   - q_as_pmf(): birth density genuinely depends on which segment a pair
#     falls in, so a mathematically correct value cannot be produced from
#     a single vector -- unlike q_ds_pmf() (below), which can.
#   - log_q_bd2_pmf()'s old (r, status) form isn't restored under that
#     name: it would have to compute birth density from source alone,
#     which is the same problem as q_as_pmf(). log_q_new_pmf() and
#     log_q_bd_pmf() keep their corrected (r_target, r_source, ...) forms
#     for the same reason and to avoid two different behaviors sharing one
#     name (was requested as "assess feasibility", not required).
#
# q_new() specifically: unlike the other five wrappers, this one does NOT
# dispatch to corrected machinery. It reproduces the exact original
# behavior -- an unconditional raw Bernoulli(N-2, p) draw, no admissibility
# check, no self-transition, no current-state argument. This is
# deliberate: the original 3-argument (lambda_r, dt, N) form has no
# current state to fall back to, and this wrapper's contract is to be the
# backward-compatible RAW generator, not a proposal kernel. The corrected
# MH sampler never calls q_new(); it calls .q_new_core(r_cur, ...)
# directly, which DOES apply the admissibility/self-transition rule.


#' @rdname proposals
#' @section Backward compatibility:
#' \code{q_new(lambda_r, dt, N)}, \code{q_bd(r_cur)}, \code{q_bd2(r_cur)},
#' \code{q_shift(r_cur)}, \code{q_as(r_cur)}, \code{q_as_pmf(r_target,
#' r_source)}, \code{q_ds(r_cur)}, and \code{q_ds_pmf(r)} restore the
#' function names (and, other than \code{q_as_pmf()}, the exact original
#' calling convention) from the pre-revision package. All but
#' \code{q_new()} dispatch to the corrected machinery above rather than
#' the original formulas; \code{q_new()} intentionally reproduces the
#' exact original (uncorrected) raw generator -- see
#' "Backward-compatibility wrappers" in the source.
#' @export
q_new <- function(lambda_r, dt, N) {
  # TRUE-original-compatible: returns the raw independent Bernoulli(N-2, p)
  # draw with NO admissibility check and NO self-transition -- exactly what
  # the pre-revision q_new(lambda_r, dt, N) did. The corrected MH sampler
  # never calls this; it calls .q_new_core(r_cur, ...) directly, which DOES
  # apply the admissibility/self-transition rule this wrapper deliberately
  # omits.
  p <- 1 - exp(-lambda_r * dt)
  stats::rbinom(N - 2, size = 1, prob = p)
}

#' @rdname proposals
#' @export
q_bd <- function(r_cur, N = length(r_cur) + 2L, K_max = NULL, min_gap = 1L) {
  .q_bd_core(r_cur, N, K_max, min_gap)
}

#' @rdname proposals
#' @export
q_bd2 <- function(r_cur, N = length(r_cur) + 2L, K_max = NULL, min_gap = 1L) {
  .q_bd2_core(r_cur, N, K_max, min_gap)
}

#' @rdname proposals
#' @export
q_as <- function(r_cur, N = length(r_cur) + 2L, K_max = NULL, min_gap = 1L) {
  .q_bd2_birth_sample(r_cur, N, K_max, min_gap)
}

#' @rdname proposals
#' @export
q_ds <- function(r_cur, N = length(r_cur) + 2L, K_max = NULL, min_gap = 1L) {
  .q_bd2_death_sample(r_cur, N, K_max, min_gap)
}

#' @rdname proposals
#' @details
#' For \code{q_as_pmf()}: both \code{r_target} and \code{r_source} are
#' required (unlike the original 1-argument \code{q_as_pmf(r)}) because
#' birth density genuinely depends on which segment a specific pair falls
#' in -- see "Backward compatibility".
#' @export
q_as_pmf <- function(r_target, r_source, N = length(r_source) + 2L,
                      K_max = NULL, min_gap = 1L) {
  exp(.q_bd2_birth_log_pmf(r_target, r_source, N, K_max, min_gap))
}

#' @rdname proposals
#' @param r (for \code{q_ds_pmf()}) the source vector. Unlike birth, every
#'   legal single-pair deletion from a given source carries the SAME
#'   density (0.5 / |D(r)|), so this can be -- and is -- computed
#'   correctly from \code{r} alone, with no target argument needed; this
#'   is the only \verb{q_*_pmf} wrapper that keeps its exact original
#'   1-argument-style signature (\code{N}/\code{K_max}/\code{min_gap} are
#'   new but default sensibly). Returns \code{0} when no legal deletion
#'   exists (fewer than 2 changepoints), matching the original's
#'   \code{if (K_r >= 2) ... else 0} guard.
#' @export
q_ds_pmf <- function(r, N = length(r) + 2L, K_max = NULL, min_gap = 1L) {
  c_r <- sum(r)
  if (c_r < 2) 0 else 0.5 / (c_r - 1L)
}
