#' Check whether a single CPLASS chain has stopped improving
#'
#' CPLASS's Metropolis-Hastings search does not converge to a stationary
#' distribution in the usual Bayesian-inference sense — it is a stochastic
#' search for the changepoint configuration that maximizes the criterion
#' \eqn{\Phi(r)}. The relevant question for a single chain is therefore not
#' "has it mixed?" but "has it stopped finding better configurations?".
#'
#' This runs \code{\link{CPLASS}} once with \code{Diagnostic = TRUE},
#' reconstructs the trace of the *accepted* state's criterion score at every
#' iteration (from `info_table`'s `score_cur`/`score_new`/`decision`
#' columns), and checks how long it has been since the running maximum last
#' improved. If nothing better has been found in the last
#' \code{window_frac} fraction of post-burn-in iterations, the chain is
#' flagged as (locally) converged.
#'
#' Important: this only tells you the *chain* has settled, not that it found
#' the *global* optimum — the criterion surface can have multiple competing
#' local maxima (see Figure 4 of the companion paper), which is exactly what
#' \code{\link{CPLASS_multistart}} is for. Use both together: this to size
#' `iter_max` correctly for a given trajectory length, and multistart to
#' guard against local-maximum disagreement between runs.
#'
#' @inheritParams CPLASS-function
#' @param window_frac fraction of post-burn-in iterations used as the "no
#'   improvement" window. Default 0.2 (last 20\%).
#' @param K_max optional maximum number of segments, passed through to the
#'   underlying \code{\link{CPLASS}} call. \code{NULL} (default) means no
#'   additional practitioner bound beyond the structural ceiling. See
#'   \code{\link{MHsearch}}.
#' @param min_gap minimum allowed segment-boundary spacing, passed through
#'   to the underlying \code{\link{CPLASS}} call. Default \code{1L}. See
#'   \code{\link{MHsearch}}.
#' @return A list with:
#'   \itemize{
#'     \item `converged`: logical.
#'     \item `best_score`: the best criterion score found.
#'     \item `last_improvement_iter`: post-burn-in iteration index at which
#'       the running maximum last improved.
#'     \item `iterations_since_improvement`, `window_size`: the check's
#'       inputs, for transparency.
#'     \item `trace`: a tibble with `iteration`, `score`, `running_max`, for
#'       plotting (see \code{\link{plot_convergence}}).
#'   }
#' @examples
#' \donttest{
#' data_file <- system.file("extdata", "Real_21_Periphery.csv", package = "cplass")
#' traj <- read.csv(data_file)
#' path <- traj[traj$index_path == 3, ]
#'
#' conv <- check_convergence(path$t, path$x, path$y, iter_max = 2000, burn_in = 200)
#' conv$converged
#' }
#' @export
check_convergence <- function(t, x, y, lambda_r = 1 / 30,
                                iter_max = 5000, burn_in = 500, s_cap = 1,
                                gamma = 1.01, speed_pen = TRUE, eta = 1, sd = NA,
                                pen = "ssic", window_frac = 0.2,
                                show_progress = FALSE, K_max = NULL, min_gap = 1L) {
  res <- CPLASS(t, x, y, lambda_r = lambda_r, iter_max = iter_max,
                burn_in = burn_in, s_cap = s_cap, gamma = gamma,
                speed_pen = speed_pen, eta = eta, Diagnostic = TRUE, sd = sd, pen = pen,
                show_progress = show_progress, K_max = K_max, min_gap = min_gap)

  if (is.null(res$pl)) {
    stop("CPLASS() failed to converge on a valid fit; cannot check convergence.",
         call. = FALSE)
  }

  it <- res$info_table
  accepted_score <- ifelse(it$decision == 1, it$score_new, it$score_cur)
  running_max <- cummax(accepted_score)
  best_score <- running_max[length(running_max)]

  last_improvement_iter <- min(which(running_max == best_score))
  n_iters <- length(running_max)
  iterations_since_improvement <- n_iters - last_improvement_iter
  window_size <- ceiling(window_frac * n_iters)
  converged <- iterations_since_improvement >= window_size

  list(
    converged = converged,
    best_score = best_score,
    last_improvement_iter = last_improvement_iter,
    iterations_since_improvement = iterations_since_improvement,
    window_size = window_size,
    n_iters = n_iters,
    trace = tibble::tibble(
      iteration = seq_len(n_iters),
      score = accepted_score,
      running_max = running_max
    ),
    cplass_result = res$pl
  )
}

#' Plot a CPLASS chain's convergence trace
#'
#' @param conv the output of \code{\link{check_convergence}}.
#' @return a ggplot object: the accepted-state criterion score per
#'   iteration (light) with the running maximum overlaid (dark), and a
#'   dashed vertical line marking the last improvement.
#' @export
plot_convergence <- function(conv) {
  verdict <- if (conv$converged) "converged (no improvement in last window)" else "NOT converged - consider raising iter_max"

  ggplot2::ggplot(conv$trace, ggplot2::aes(x = .data$iteration)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$score), col = "gray70", linewidth = 0.3) +
    ggplot2::geom_step(ggplot2::aes(y = .data$running_max), col = "steelblue", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = conv$last_improvement_iter, linetype = "dashed", col = "orange") +
    ggplot2::labs(
      title = paste0("CPLASS chain convergence: ", verdict),
      subtitle = paste0(
        "best score = ", round(conv$best_score, 2),
        " | last improved at iter ", conv$last_improvement_iter,
        " of ", conv$n_iters
      ),
      x = "post-burn-in iteration", y = expression(Phi(r))
    ) +
    ggplot2::theme_classic()
}
