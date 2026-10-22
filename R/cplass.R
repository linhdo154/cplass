#' Run CPLASS on a single 2D trajectory
#'
#' Fits the CPLASS model (Section 2 of the companion paper) to a single
#' particle trajectory given by observed times \code{t} and coordinates
#' \code{x}, \code{y}: runs the Metropolis-Hastings search
#' (\code{\link{MHsearch}}) over changepoint configurations, selects the
#' configuration with the highest criterion score, and returns the
#' corresponding continuous piecewise-linear segmentation.
#'
#' Same signature and return value as the original research-script
#' \code{CPLASS()}. See \code{\link{MHsearch}} for details on what changed
#' internally for speed.
#'
#' @param t numeric vector of observation times.
#' @param x numeric vector of x-coordinates.
#' @param y numeric vector of y-coordinates.
#' @param lambda_r rate parameter for the Type-1 (independent) changepoint
#'   proposal. Default \code{1/30}.
#' @param iter_max number of MCMC iterations. Default 5000.
#' @param burn_in number of initial iterations discarded. Default 500.
#' @param s_cap maximum segment speed with no penalty. Default 1.
#' @param gamma exponent for the sSIC penalty. Default 1.01 (see Section
#'   3.1.1 of the companion paper for the rationale).
#' @param speed_pen logical; whether to activate the speed penalty.
#'   Default \code{TRUE}.
#' @param eta non-negative, dimensionless strength of the relative-excess
#'   speed penalty. Default 1. Setting \code{eta = 0} recovers the criterion
#'   without speed regularization.
#' @param Diagnostic logical; if \code{TRUE}, also return the MCMC
#'   diagnostic trace (`info_table`, `update_info`). Default \code{FALSE}.
#' @param sd optional known noise sd.
#' @param pen one of \code{"ssic"}, \code{"aicc"}, \code{"hybrid"}.
#' @param show_progress logical; display a text progress bar. Default
#'   \code{TRUE}.
#' @param patience if not \code{NULL} (default), stop the search early once
#'   the criterion score hasn't improved for this many consecutive
#'   iterations, rather than always running the full `iter_max`. See
#'   \code{\link{MHsearch}}'s "Early stopping" section for details and
#'   caveats (in particular: this doesn't protect against local maxima —
#'   pair it with \code{\link{CPLASS_multistart}} if that's a concern).
#' @param K_max optional maximum number of segments (the manuscript's
#'   practitioner-specified \eqn{\bar k} in Eq. 2.6). \code{NULL} (default)
#'   means no additional practitioner bound is imposed beyond the
#'   structural \eqn{K(r) \le n-2} ceiling -- the manuscript intentionally
#'   leaves \eqn{\bar k} to the practitioner. See \code{\link{MHsearch}}.
#' @param min_gap minimum allowed segment-boundary spacing, in
#'   observation-index units. Default \code{1L} (structurally non-binding).
#'   See \code{\link{MHsearch}}.
#' @return A list with `segments_inferred` and `path_inferred` tibbles (see
#'   \code{\link{piecewise_linear_con}}, \code{old_version = FALSE}), plus
#'   `best_score`: the criterion value \eqn{\Phi(r)} of the selected
#'   configuration, useful for comparing chains (see
#'   \code{\link{CPLASS_multistart}}, \code{\link{check_convergence}}), and
#'   `stopped_early` / `iterations_used` (only meaningful if `patience` was
#'   set). If \code{Diagnostic = TRUE}, also includes `info_table` and
#'   `update_info`. Returns \code{pl = NULL} (wrapped the same way) if the
#'   sampler fails to converge on a valid fit after retries.
#' @examples
#' \donttest{
#' data_file <- system.file("extdata", "Real_21_Periphery.csv", package = "cplass")
#' traj <- read.csv(data_file)
#' path <- traj[traj$index_path == 3, ]
#'
#' fit <- CPLASS(path$t, path$x, path$y, iter_max = 2000, burn_in = 200,
#'                show_progress = FALSE)
#' fit$segments_inferred
#' }
#' @name CPLASS-function
#' @export
CPLASS <- function(t, x, y, lambda_r = 1 / 30,
                    iter_max = 5000, burn_in = 500, s_cap = 1,
                    gamma = 1.01, speed_pen = TRUE, eta = 1,
                    Diagnostic = FALSE,
                    sd = NA, pen = "ssic", show_progress = TRUE,
                    patience = NULL, K_max = NULL, min_gap = 1L) {
  dt <- min(diff(t))
  speed_control <- if (isTRUE(speed_pen)) 1 else 0

  max_retries <- 20
  attempt <- 1
  pl <- list()
  info_table <- NULL
  update_info <- NULL

  repeat {
    MCMC <- tryCatch(
      MHsearch(t, x, y, dt, lambda_r = lambda_r,
               iter_max = iter_max, burn_in = burn_in,
               s_cap = s_cap, gamma = gamma,
               speed_control = speed_control, eta = eta, sd = sd, pen = pen,
               show_progress = show_progress, patience = patience,
               K_max = K_max, min_gap = min_gap),
      error = function(e) e
    )

    if (!inherits(MCMC, "error")) {
      Final_MCMC <- FinalMH(MCMC$cps_list)
      index.max <- which.max(Final_MCMC$criterion_score)
      pl <- piecewise_linear_con(t, x, y, Final_MCMC$this_cp[[index.max]],
                                   old_version = FALSE)
      pl$best_score <- Final_MCMC$criterion_score[index.max]
      pl$stopped_early <- MCMC$stopped_early
      pl$iterations_used <- MCMC$iterations_used
      if (isTRUE(Diagnostic)) {
        info_table <- MCMC$info_table
        update_info <- MCMC$update_info
      }
      break
    } else {
      message("MHsearch failed, retrying... attempt ", attempt)
      attempt <- attempt + 1
      if (attempt > max_retries) {
        message("MHsearch failed too many times. Returning pl = NULL.")
        pl <- NULL
        break
      }
    }
  }

  if (isTRUE(Diagnostic)) {
    list(pl = pl, info_table = info_table, update_info = update_info)
  } else {
    pl
  }
}
