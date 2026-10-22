#' Metropolis-Hastings search over the changepoint space
#'
#' Runs the tailored MH sampler described in Section 2.4 of the companion
#' paper for a single 2D trajectory. This has the same signature and return
#' structure as the original research-script \code{MHsearch()}. See
#' "Performance" and "Proposal mechanism" below for how the internals
#' differ from that original script.
#'
#' @section Performance:
#' The original implementation recomputed the full piecewise-linear fit
#' (\code{\link{piecewise_linear_con}}) up to three times per iteration:
#' once for the current state, once for the proposal, and once more after
#' the accept/reject step to populate the trace. Since the fit is a pure
#' function of \code{(t, x, y, changepoints)}, the current state's fit
#' never needs to be recomputed (it's identical to the state fit at the end
#' of the previous iteration), and the post-decision refit is always
#' identical to either the already-computed current or proposal fit. This
#' version fits the proposal once per iteration and carries the accepted
#' fit forward, cutting the number of `piecewise_linear_con()`/`CS()` calls
#' from ~3 to 1 per iteration. The diagnostic trace (`info_table`) is also
#' built from preallocated vectors instead of row-by-row `bind_rows()`,
#' which previously made the trace step scale quadratically in
#' `iter_max`.
#'
#' @section Proposal mechanism:
#' The changepoint proposal mechanism (\code{\link{proposals}}) mixes four
#' move types with fixed weights (1/4, 1/8, 1/8, 1/2) at every state, and
#' evaluates the forward and reverse proposal densities for the same move
#' type that generated the proposal. This is a deliberate correction, not
#' a performance refactor: an earlier version of this mechanism
#' redistributed the move-type weights at boundary/saturated states and
#' could re-derive a different move type for the reverse-direction density
#' than the one that actually generated the proposal, which does not
#' satisfy detailed balance. Results are therefore not expected to match
#' that earlier mechanism's accept/reject decisions, even given the same
#' random-draw stream.
#'
#' @section Early stopping:
#' If \code{patience} is set, \code{iter_max} becomes a safety ceiling
#' rather than a fixed target: once the running-maximum criterion score
#' (post-burn-in) hasn't improved for \code{patience} consecutive
#' iterations, the chain stops early. This is the same "no improvement"
#' idea as \code{\link{check_convergence}}, applied live inside the loop so
#' you don't pay for iterations a path doesn't need. It does not protect
#' against local maxima in the criterion surface (see
#' \code{\link{CPLASS_multistart}} for that) — a chain can plateau on a
#' local max just as easily as a global one.
#'
#' @inheritParams piecewise_linear_con
#' @param dt observation time step.
#' @param lambda_r rate parameter for the Type-1 (independent) changepoint
#'   proposal.
#' @param iter_max number of MCMC iterations (an upper bound if `patience`
#'   is set; see "Early stopping").
#' @param burn_in number of initial iterations discarded from the trace.
#' @param s_cap maximum segment speed with no penalty.
#' @param eta non-negative, dimensionless strength of the relative-excess
#'   speed penalty. Default 1; \code{eta = 0} removes this regularization.
#' @param gamma exponent for the sSIC penalty.
#' @param speed_control 1 to activate the speed penalty, 0 to deactivate.
#' @param sd optional known noise sd (see \code{\link{loglikelihood}}).
#' @param pen one of \code{"ssic"}, \code{"aicc"}, \code{"hybrid"}.
#' @param show_progress logical; display a text progress bar (default
#'   \code{TRUE}, matching the original behavior).
#' @param patience if not \code{NULL} (default), stop early once the
#'   post-burn-in running-maximum criterion score hasn't improved for this
#'   many consecutive iterations. See "Early stopping".
#' @param K_max optional maximum number of segments (the manuscript's
#'   practitioner-specified \eqn{\bar k} in Eq. 2.6); passed through to
#'   \code{\link{CS}} and the proposal mechanism (\code{\link{proposals}}).
#'   \code{NULL} (default) means no additional practitioner bound is
#'   imposed beyond the structural \eqn{K(r) \le n-2} ceiling -- the
#'   manuscript intentionally leaves \eqn{\bar k} to the practitioner.
#' @param min_gap minimum allowed segment-boundary spacing (observation-
#'   index units), passed through to \code{\link{CS}} and the proposal
#'   mechanism. Default \code{1L} (structurally non-binding).
#' @return A list with `cps_list` (post-burn-in trace), `info_table`
#'   (post-burn-in diagnostics), `update_info` (acceptance rate), and
#'   `stopped_early` / `iterations_used` (whether/when patience triggered;
#'   `iterations_used` equals `iter_max` if it never did).
#' @export
MHsearch <- function(t, x, y, dt, lambda_r = 1 / 30,
                      iter_max = 5000, burn_in = 500, s_cap = 1,
                      gamma = 1.01, speed_control = 0, eta = 1, sd = NA,
                      pen = "ssic",
                      show_progress = TRUE, patience = NULL,
                      K_max = NULL, min_gap = 1L) {
  .validate_K_max_min_gap(K_max, min_gap)
  N <- length(t)
  penalty_coef <- max(4, log(N))^gamma

  fit_state <- function(r) {
    CS(t, x, y, r, penalty_coef, s_cap, speed_control, eta = eta, sd = sd,
       pen = pen, gamma = gamma, K_max = K_max, min_gap = min_gap)
  }

  r0 <- .q_new_core(integer(N - 2L), N, lambda_r, dt, K_max, min_gap)
  cur <- fit_state(r0)
  attempts <- 0L
  while (!isTRUE(cur$logical) && attempts < 100L) {
    r0 <- .q_new_core(integer(N - 2L), N, lambda_r, dt, K_max, min_gap)
    cur <- fit_state(r0)
    attempts <- attempts + 1L
  }
  if (!isTRUE(cur$logical)) {
    stop("Could not find a valid initial changepoint configuration after ",
         "100 attempts.", call. = FALSE)
  }

  if (burn_in >= iter_max) {
    stop(
      "`burn_in` (", burn_in, ") must be smaller than `iter_max` (", iter_max,
      "); the original script silently produced an empty post-burn-in trace ",
      "in this case, which crashed downstream in FinalMH() with a cryptic ",
      "'subscript out of bounds' error. Increase `iter_max` or lower `burn_in`.",
      call. = FALSE
    )
  }

  cps_list <- vector("list", iter_max)

  llh_cur <- llh_new <- score_cur <- score_new <- ncp_cur <- ncp_new <-
    logA_v <- pen_cur <- pen_new <- est_sigma_cur <- est_sigma_new <-
    decision_v <- numeric(iter_max)

  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = iter_max, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  running_best <- -Inf
  last_improve_at <- burn_in
  actual_iters <- iter_max
  stopped_early <- FALSE

  for (count in seq_len(iter_max)) {
    pp <- proposal_function(r0, N, lambda_r, dt, K_max = K_max, min_gap = min_gap)
    r_prop <- pp$r_prop
    move_type <- pp$move_type

    if (identical(r_prop, r0)) {
      # Self-transition: the proposal mechanism already determined no
      # legal off-diagonal move of this type is available from the current
      # state (or an admissible draw happened to equal it). No refit, no
      # proposal-density evaluation, and in particular no Type-1 self-loop
      # DP (.bernoulli_admissible_mass()) is invoked here -- that exact
      # diagonal density is only ever computed when something explicitly
      # queries pproposal()/log_q_new_pmf() with r_target == r_source
      # (the unit tests and exhaustive detailed-balance checks do this;
      # this hot loop never does).
      r1 <- r0
      new_state <- cur
      prop <- cur
      logA <- 0
    } else {
      prop <- fit_state(r_prop)

      if (!isTRUE(prop$logical) || !isTRUE(cur$logical)) {
        # Hard rejection: an invalid fit (e.g. rank-deficient) can never be
        # accepted, independent of its criterion score. A criterion score
        # that happens to equal exactly 0 is NOT a reason to reject on its
        # own (exp(0) = 1 is a perfectly valid contribution to the MH
        # ratio) -- only `logical = FALSE` identifies an invalid state.
        logA <- -Inf
        accept <- FALSE
      } else {
        # CRITICAL: forward and reverse proposal densities are evaluated for
        # the SAME move_type that actually generated r_prop -- never
        # re-derived from a fresh u draw or from which state is "given".
        log_q_fwd <- pproposal(move_type, r_prop, r0, N, lambda_r, dt, K_max, min_gap)
        log_q_rev <- pproposal(move_type, r0, r_prop, N, lambda_r, dt, K_max, min_gap)
        logA <- (prop$s - cur$s) + (log_q_rev - log_q_fwd)
        accept <- is.finite(logA) && (logA >= 0 || stats::runif(1) < exp(logA))
      }

      if (accept) {
        r1 <- r_prop
        new_state <- prop
      } else {
        r1 <- r0
        new_state <- cur
      }
    }

    pla1 <- new_state$pla
    cps1 <- pla1$this_cp

    cps_list[[count]] <- list(
      r = r1,
      u = pla1$u_x,
      v = pla1$v_y,
      this_cp = cps1,
      seg_vel = pla1$path_segvel,
      seg_time = pla1$path_segtime,
      eta = pla1$path_segeta,
      criterion_score = new_state$s,
      RSS = pla1$path_RSS,
      r_prop = r_prop
    )

    llh_cur[count] <- cur$llh
    llh_new[count] <- if (isTRUE(prop$logical)) prop$llh else 0
    score_cur[count] <- cur$s
    score_new[count] <- if (isTRUE(prop$logical)) prop$s else 0
    ncp_cur[count] <- sum(r0)
    ncp_new[count] <- sum(r_prop)
    logA_v[count] <- if (is.finite(logA)) logA else 0
    pen_cur[count] <- cur$p
    pen_new[count] <- if (isTRUE(prop$logical)) prop$p else 0
    est_sigma_cur[count] <- sqrt(cur$pla$path_RSS / (2 * N))
    est_sigma_new[count] <- if (isTRUE(prop$logical)) sqrt(prop$pla$path_RSS / (2 * N)) else 0
    decision_v[count] <- as.numeric(any(r1 != r0))

    r0 <- r1
    cur <- new_state

    if (show_progress) utils::setTxtProgressBar(pb, count)

    if (!is.null(patience) && count > burn_in) {
      if (new_state$s > running_best) {
        running_best <- new_state$s
        last_improve_at <- count
      }
      if ((count - last_improve_at) >= patience) {
        actual_iters <- count
        stopped_early <- TRUE
        break
      }
    }
  }

  # Truncate to the iterations actually run (no-op if we ran the full iter_max)
  cps_list <- cps_list[seq_len(actual_iters)]
  llh_cur <- llh_cur[seq_len(actual_iters)]
  llh_new <- llh_new[seq_len(actual_iters)]
  score_cur <- score_cur[seq_len(actual_iters)]
  score_new <- score_new[seq_len(actual_iters)]
  ncp_cur <- ncp_cur[seq_len(actual_iters)]
  ncp_new <- ncp_new[seq_len(actual_iters)]
  logA_v <- logA_v[seq_len(actual_iters)]
  pen_cur <- pen_cur[seq_len(actual_iters)]
  pen_new <- pen_new[seq_len(actual_iters)]
  est_sigma_cur <- est_sigma_cur[seq_len(actual_iters)]
  est_sigma_new <- est_sigma_new[seq_len(actual_iters)]
  decision_v <- decision_v[seq_len(actual_iters)]

  info_table <- tibble::tibble(
    llh_cur = llh_cur,
    llh_new = llh_new,
    score_cur = score_cur,
    score_new = score_new,
    ncp_cur = ncp_cur,
    ncp_new = ncp_new,
    logA = logA_v,
    pen_cur = pen_cur,
    pen_new = pen_new,
    est_sigma_cur = est_sigma_cur,
    est_sigma_new = est_sigma_new,
    ind_iteration = seq_len(actual_iters),
    decision = decision_v
  )

  keep <- seq.int(burn_in + 1L, actual_iters)
  it <- info_table[keep, ]

  update_info <- tibble::tibble(
    accept_rate = sum(it$decision) / length(it$decision),
    reject_rate = 1 - sum(it$decision) / length(it$decision)
  )

  list(
    cps_list = cps_list[keep],
    info_table = it,
    update_info = update_info,
    stopped_early = stopped_early,
    iterations_used = actual_iters
  )
}

#' Collect a Metropolis-Hastings trace into flat vectors/lists
#'
#' Drop-in, vectorized replacement for the original \code{FinalMH()}. The
#' original also computed a `loglikelihood` vector that was never included
#' in its return value (dead code); this version omits that computation.
#'
#' @param bcp_new the `cps_list` element of \code{\link{MHsearch}}'s output.
#' @return A list of parallel vectors/lists extracted from the trace:
#'   `r`, `u`, `v`, `this_cp`, `seg_vel`, `seg_time`, `eta`,
#'   `criterion_score`, `RSS`, `num_cp`, `r_prop`.
#' @export
FinalMH <- function(bcp_new) {
  list(
    r = lapply(bcp_new, `[[`, "r"),
    u = lapply(bcp_new, `[[`, "u"),
    v = lapply(bcp_new, `[[`, "v"),
    this_cp = lapply(bcp_new, `[[`, "this_cp"),
    seg_vel = lapply(bcp_new, `[[`, "seg_vel"),
    seg_time = lapply(bcp_new, `[[`, "seg_time"),
    eta = vapply(bcp_new, `[[`, numeric(1), "eta"),
    criterion_score = vapply(bcp_new, `[[`, numeric(1), "criterion_score"),
    RSS = vapply(bcp_new, `[[`, numeric(1), "RSS"),
    num_cp = vapply(bcp_new, function(z) length(z$this_cp), integer(1)),
    r_prop = lapply(bcp_new, `[[`, "r_prop")
  )
}
