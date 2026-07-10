#' Pathwise mean squared displacement (MSD)
#'
#' @param path a path object with `path$path$x`, `path$path$y`,
#'   `path$path$t`.
#' @param use_pct fraction of the trajectory length to use as the maximum
#'   lag. Default 0.2.
#' @return a tibble with columns `t_lag` and `msd`.
#' @export
pathwise_MSD <- function(path, use_pct = 0.2) {
  x <- path$path$x
  y <- path$path$y
  t <- path$path$t

  n <- length(t)
  n_lag <- floor(use_pct * n)
  t_lag <- t[seq_len(n_lag)]

  # vectorized: for lag `tt`, average squared displacement over all pairs
  # (i, i+tt) at once, instead of a nested for loop over lag x offset.
  msd <- vapply(seq_len(n_lag), function(tt) {
    idx <- seq_len(n - tt)
    mean((x[idx + tt] - x[idx])^2 + (y[idx + tt] - y[idx])^2) / 4
  }, numeric(1))

  tibble::tibble(t_lag = t_lag, msd = msd)
}

#' Ensemble (group-averaged) MSD
#'
#' @param msd_df a tibble stacked across paths (e.g. via
#'   \code{dplyr::bind_rows(msd_list, .id = "id")}), with columns `t_lag`
#'   and `msd`.
#' @return a tibble with columns `t_lag` and `e_msd`.
#' @details
#' NOTE: the original script's `ensemble_MSD()` grouped by a column named
#' `t`, but was always called on plain lists of `pathwise_MSD()` outputs
#' (which use `t_lag`, not `t`/`id`) rather than on a combined data frame,
#' so `plot_msd()`'s calls to `ensemble_MSD(msd_list)` would have errored
#' as originally written. This version is made internally consistent
#' (grouping by `t_lag` on a properly stacked data frame), but the exact
#' intended grouping/labeling should be confirmed.
#' @export
ensemble_MSD <- function(msd_df) {
  msd_df |>
    dplyr::group_by(.data$t_lag) |>
    dplyr::summarise(e_msd = mean(.data$msd), .groups = "drop")
}

#' Plot pathwise and ensemble MSD curves for a group, split by activity
#'
#' @param msd_list a list of per-path MSD tibbles (`pathwise_MSD()` output).
#' @param is_active logical vector, same length as `msd_list`.
#' @param group_name label used in the plot title.
#' @return invisibly, the overall ensemble MSD tibble.
#' @export
plot_msd <- function(msd_list, is_active, group_name = "") {
  emsd <- ensemble_MSD(dplyr::bind_rows(msd_list, .id = "id"))
  full_msd_t <- emsd$t_lag
  dt <- full_msd_t[2] - full_msd_t[1]
  active_list <- which(is_active)
  not_active <- which(!is_active)

  plot(1,
    type = "l", xlim = c(dt, max(full_msd_t)), ylim = c(0.001, 1),
    log = "xy", xlab = "Time (s)", ylab = "Microns^2",
    main = paste0("Group MSD: ", group_name)
  )

  for (i in not_active) {
    graphics::lines(msd_list[[i]]$t_lag, msd_list[[i]]$msd, type = "l", col = "pink")
  }
  for (i in active_list) {
    graphics::lines(msd_list[[i]]$t_lag, msd_list[[i]]$msd, type = "l", col = "darkgray")
  }

  emsd_not_active <- ensemble_MSD(dplyr::bind_rows(msd_list[not_active], .id = "id"))
  emsd_active <- ensemble_MSD(dplyr::bind_rows(msd_list[active_list], .id = "id"))

  graphics::lines(emsd_not_active$t_lag, emsd_not_active$e_msd, lwd = 2, col = "red")
  graphics::lines(emsd_active$t_lag, emsd_active$e_msd, lwd = 2, col = "black")
  graphics::lines(emsd$t_lag, emsd$e_msd, lwd = 2, col = "blue")

  graphics::legend("topleft",
    legend = c(
      "Ensemble MSD", "Pathwise MSD, Not Active", "Ensemble MSD, Not Active",
      "Pathwise MSD, Active", "Ensemble MSD, Active"
    ),
    col = c("blue", "pink", "red", "gray", "black"),
    lwd = c(2, 1, 2, 1, 2)
  )

  invisible(emsd)
}
