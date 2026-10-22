#' Run CPLASS from multiple independent random starts and keep the best
#'
#' Runs \code{\link{CPLASS}} \code{n_starts} times with different seeds and
#' returns the run with the highest criterion score
#' (\code{best_score}), along with a summary of how much the runs agreed —
#' the multi-chain analogue of \code{\link{check_convergence}}'s
#' single-chain check. Large disagreement in \code{all_scores} or
#' \code{all_n_segments} indicates the criterion surface has multiple
#' comparably-good local maxima for this trajectory (see Figure 4 of the
#' companion paper) — in that case, the *Cumulative Speed Allocation*
#' summary (\code{\link{compute_csa}}) is a more robust downstream quantity
#' than any single run's segment count.
#'
#' @section Patience + multistart, together:
#' By default (\code{patience = "auto"}), each of the \code{n_starts}
#' chains uses \code{\link{MHsearch}}'s patience-based early stopping
#' (ceiling(0.3 * iter_max) iterations without improvement), so no single
#' chain pays for iterations it doesn't need. Multistart's job is to make
#' that safe: because each chain has an independent random trajectory
#' through the search space, it's very unlikely that two different chains
#' plateau at exactly the same (possibly suboptimal) configuration, so
#' taking the max across chains recovers what any one early-stopped chain
#' might have missed. Empirically (see the package's benchmark scripts),
#' `n_starts = 5` with auto patience matches the best score of 5 full-length
#' (no-patience) runs almost exactly, at a fraction of the time. Set
#' \code{patience = NULL} to disable early stopping entirely and run every
#' chain to the full `iter_max` (slower, and in our tests did not find
#' better results than the default — but available if you want to verify
#' that yourself on your own data).
#'
#' @inheritParams CPLASS-function
#' @param n_starts number of independent chains to run. Default 5.
#' @param patience early-stopping patience passed to each chain's
#'   \code{\link{MHsearch}} call. Default \code{"auto"}: uses
#'   \code{ceiling(0.3 * iter_max)}. Set to \code{NULL} to disable early
#'   stopping (every chain runs the full `iter_max`), or supply a number to
#'   use that patience directly. See "Patience + multistart, together".
#' @param seeds optional integer vector of length `n_starts`; if `NULL`
#'   (default), seeds are drawn from the current RNG state, and are
#'   returned in the output so the winning run can be reproduced exactly
#'   with \code{set.seed(result$best_seed)}.
#' @param n_cores number of cores to use via \code{parallel::mclapply()}
#'   (not available on Windows; falls back to sequential with a message).
#'   Default 1 (sequential).
#' @param K_max optional maximum number of segments, passed through to
#'   every chain's \code{\link{CPLASS}} call. \code{NULL} (default) means
#'   no additional practitioner bound beyond the structural ceiling. See
#'   \code{\link{MHsearch}}.
#' @param min_gap minimum allowed segment-boundary spacing, passed through
#'   to every chain's \code{\link{CPLASS}} call. Default \code{1L}. See
#'   \code{\link{MHsearch}}.
#' @return A list with:
#'   \itemize{
#'     \item `best`: the winning \code{\link{CPLASS}} output (highest
#'       `best_score`).
#'     \item `best_seed`: the seed that produced it — pass to
#'       `set.seed()` before re-running `CPLASS()` to reproduce it exactly.
#'     \item `all_scores`, `all_n_segments`, `seeds`: per-chain results, for
#'       inspecting agreement across starts.
#'     \item `score_spread`: `max(all_scores) - min(all_scores)`; small
#'       values indicate the starts agree on (approximately) the same
#'       optimum.
#'     \item `patience_used`: the actual patience value applied to every
#'       chain (resolved from `"auto"` if applicable), for transparency.
#'   }
#' @examples
#' \donttest{
#' data_file <- system.file("extdata", "Real_21_Periphery.csv", package = "cplass")
#' traj <- read.csv(data_file)
#' path <- traj[traj$index_path == 3, ]
#'
#' ms <- CPLASS_multistart(path$t, path$x, path$y, iter_max = 2000, burn_in = 200,
#'                          n_starts = 3, show_progress = FALSE)
#' ms$all_scores
#' ms$score_spread
#' }
#' @export
CPLASS_multistart <- function(t, x, y, lambda_r = 1 / 30,
                                iter_max = 5000, burn_in = 500, s_cap = 1,
                                gamma = 1.01, speed_pen = TRUE, eta = 1, sd = NA,
                                pen = "ssic", n_starts = 5, patience = "auto",
                                seeds = NULL, n_cores = 1, show_progress = FALSE,
                                K_max = NULL, min_gap = 1L) {
  if (identical(patience, "auto")) {
    patience_used <- ceiling(0.3 * iter_max)
  } else {
    patience_used <- patience
  }

  if (is.null(seeds)) {
    seeds <- sample.int(.Machine$integer.max / 2, n_starts)
  }
  if (length(seeds) != n_starts) {
    stop("`seeds` must have length `n_starts` (", n_starts, ").", call. = FALSE)
  }

  run_one <- function(s) {
    set.seed(s)
    CPLASS(t, x, y, lambda_r = lambda_r, iter_max = iter_max, burn_in = burn_in,
           s_cap = s_cap, gamma = gamma, speed_pen = speed_pen, eta = eta, sd = sd,
           pen = pen, show_progress = show_progress, patience = patience_used,
           K_max = K_max, min_gap = min_gap)
  }

  if (n_cores > 1) {
    if (.Platform$OS.type == "windows") {
      message("n_cores > 1 is not supported on Windows via parallel::mclapply(); running sequentially.")
      results <- lapply(seeds, run_one)
    } else {
      results <- parallel::mclapply(seeds, run_one, mc.cores = n_cores)
    }
  } else {
    results <- lapply(seeds, run_one)
  }

  failed <- vapply(results, is.null, logical(1))
  if (all(failed)) {
    stop("All ", n_starts, " starts failed to converge on a valid fit.", call. = FALSE)
  }
  if (any(failed)) {
    message(sum(failed), " of ", n_starts, " starts failed to converge and were dropped.")
    results <- results[!failed]
    seeds <- seeds[!failed]
  }

  scores <- vapply(results, function(r) r$best_score, numeric(1))
  n_segments <- vapply(results, function(r) nrow(r$segments_inferred), integer(1))
  best_idx <- which.max(scores)

  list(
    best = results[[best_idx]],
    best_seed = seeds[best_idx],
    all_scores = scores,
    all_n_segments = n_segments,
    seeds = seeds,
    score_spread = max(scores) - min(scores),
    patience_used = patience_used
  )
}
