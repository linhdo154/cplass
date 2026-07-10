#' Cumulative Speed Allocation (CSA)
#'
#' For every speed in \code{speed_mesh}, computes the duration-weighted
#' proportion of time spent at speeds less than or equal to that value
#' (Cook et al. 2025; Section 3.1.2 of the companion paper).
#'
#' @param segments_summary a data frame/tibble with columns `durations` and
#'   `speeds`.
#' @param speed_mesh numeric vector of speed values at which to evaluate
#'   the CSA.
#' @return a tibble with columns `s` (the mesh) and `csa` (cumulative
#'   proportion of time).
#' @examples
#' segs <- data.frame(durations = c(2, 1, 3), speeds = c(0.01, 0.5, 0.05))
#' compute_csa(segs, speed_mesh = seq(0, 0.6, by = 0.1))
#' @export
compute_csa <- function(segments_summary, speed_mesh) {
  total_time <- sum(segments_summary$durations)
  # vectorized over the mesh: for each s in speed_mesh, sum durations of
  # segments slower than s, normalized by total time
  wecdf <- vapply(speed_mesh, function(s) {
    sum(segments_summary$durations * (segments_summary$speeds < s)) / total_time
  }, numeric(1))

  tibble::tibble(s = speed_mesh, csa = wecdf)
}

#' Theoretical CSA under a two-state (stationary/motile) kinetic model
#'
#' @param theta a list with elements `speed_alpha`, `speed_beta`,
#'   `dist_avg`, `p_MM`, `dur_stationary_avg`, `DEPENDENT_SPEED_DUR`.
#' @param s_mesh numeric vector of speeds.
#' @return numeric vector, the theoretical CSA evaluated at `s_mesh`.
#' @export
csa_theoretical <- function(theta, s_mesh) {
  alpha <- theta$speed_alpha
  beta <- theta$speed_beta
  c_rate <- 1 / theta$dist_avg
  num_M_avg <- 1 / (1 - theta$p_MM)
  dur_S_avg <- theta$dur_stationary_avg

  if (theta$DEPENDENT_SPEED_DUR) {
    dur_M_avg <- beta / (c_rate * (alpha - 1))
    exp_delta_I <- dur_S_avg + num_M_avg * dur_M_avg * stats::pgamma(s_mesh * 1000, shape = alpha - 1, rate = beta)
  } else {
    dur_M_avg <- theta$dist_avg / (alpha / beta)
    exp_delta_I <- dur_S_avg + num_M_avg * dur_M_avg * stats::pgamma(s_mesh * 1000, shape = alpha, rate = beta)
  }

  exp_delta_T <- dur_S_avg + num_M_avg * dur_M_avg
  exp_delta_I / exp_delta_T
}

#' @keywords internal
#' @noRd
.summarize_segments_impl <- function(path_list, group_label, field, state_col = "states") {
  out <- vector("list", length(path_list))
  for (i in seq_along(path_list)) {
    seg <- path_list[[i]][[field]][, c("durations", state_col, "speeds")]
    names(seg)[names(seg) == state_col] <- "states"
    seg$label <- group_label
    seg$path_id <- i
    out[[i]] <- seg
  }
  dplyr::bind_rows(out)
}

#' Stack segment summaries from a list of true (ground-truth) paths
#' @param path_list a list of path objects, each with a `$segments` tibble.
#' @param group_label label to attach to every row (e.g. experimental group).
#' @return a tibble with columns `durations`, `states`, `speeds`, `label`,
#'   `path_id`.
#' @export
summarize_segments <- function(path_list, group_label) {
  .summarize_segments_impl(path_list, group_label, field = "segments")
}

#' Stack segment summaries from a list of CPLASS-inferred paths
#' @inheritParams summarize_segments
#' @export
summarize_segments_inferred <- function(path_list, group_label) {
  .summarize_segments_impl(path_list, group_label, field = "segments_inferred")
}

#' Stack segment summaries using a threshold-based test-state column
#' @inheritParams summarize_segments
#' @export
summarize_segments_test <- function(path_list, group_label) {
  .summarize_segments_impl(path_list, group_label,
                            field = "segments_inferred", state_col = "states_test")
}

#' Stack segment summaries from inferred paths (cutoff-labeled states)
#' @inheritParams summarize_segments
#' @export
summarize_segments_cutoff <- function(path_list, group_label) {
  .summarize_segments_impl(path_list, group_label, field = "segments_inferred")
}
