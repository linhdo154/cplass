col_actual <- "orange"
col_inferred <- "steelblue"
col_stationary <- "red3"
col_motile <- "green4"
col_stationary_inferred <- "pink"
col_motile_inferred <- "green"

#' @keywords internal
#' @noRd
.title_index <- function(path_index) {
  # The original scripts referenced a bare `i` in plot titles that was
  # never a formal argument, relying instead on whatever `i` happened to
  # be left over in the calling environment (e.g. a `for (i in ...)` loop
  # around the call). That works but is fragile outside that exact usage
  # pattern. This keeps the same default behavior (falls back to a caller's
  # `i`, or "" if none exists) while also accepting an explicit
  # `path_index` argument as a safer alternative.
  if (!is.null(path_index)) {
    return(path_index)
  }
  get0("i", envir = parent.frame(2), ifnotfound = "")
}

#' Compute plot bounds for a single inferred path
#' @param path_info a CPLASS output list (with `path_inferred`).
#' @export
get_one_plot_size <- function(path_info) {
  path <- path_info$path_inferred
  dist_max <- max(max(path$x) - min(path$x), max(path$y) - min(path$y))
  t_max <- max(path$t) - min(path$t)
  tibble::tibble(xy_width = 1.2 * dist_max, t_begin = 0, t_end = t_max)
}

#' Compute shared plot bounds across a list of ground-truth paths
#' @param path_list a list of path objects, each with `$path`.
#' @export
get_plot_size <- function(path_list) {
  dist_max <- vapply(path_list, function(p) {
    path <- p$path
    max(max(path$x) - min(path$x), max(path$y) - min(path$y))
  }, numeric(1))
  t_max <- vapply(path_list, function(p) {
    path <- p$path
    max(path$t) - min(path$t)
  }, numeric(1))
  tibble::tibble(xy_width = 1.2 * max(dist_max), t_begin = 0, t_end = max(t_max))
}

#' Compute shared plot bounds across a list of inferred paths
#' @param path_list a list of path objects, each with `$path_inferred`.
#' @export
get_plot_size_inferred <- function(path_list) {
  dist_max <- vapply(path_list, function(p) {
    path <- p$path_inferred
    max(max(path$x) - min(path$x), max(path$y) - min(path$y))
  }, numeric(1))
  t_max <- vapply(path_list, function(p) {
    path <- p$path_inferred
    max(path$t) - min(path$t)
  }, numeric(1))
  tibble::tibble(xy_width = 1.2 * max(dist_max), t_begin = 0, t_end = max(t_max))
}

#' Four-panel dashboard for a single CPLASS-segmented trajectory
#'
#' Shows the x-y path with inferred anchor overlay, the x(t) and y(t) time
#' series with motile/stationary shading, and a duration-vs-speed scatter
#' of inferred segments.
#'
#' @param path_info CPLASS output for one path (`path_inferred`,
#'   `segments_inferred`).
#' @param xy_width width of the x-y panel; computed automatically if `NA`.
#' @param t_lim time-axis limits; computed automatically if `NA`.
#' @param motor label used in the plot title.
#' @param max_speed upper limit for the speed axis in the duration/speed panel.
#' @param state_shaded logical; shade motile/stationary segments. Default `TRUE`.
#' @param title_ind logical; include a path index in the title. Default `TRUE`.
#' @param show_time_changes logical; draw vertical lines at changepoints
#'   instead of state shading. Default `FALSE`.
#' @param path_index optional explicit path index/label for the title (see
#'   Details in package overview about the original's implicit `i` lookup).
#' @return a `patchwork` composed ggplot object.
#' @examples
#' \donttest{
#' data_file <- system.file("extdata", "Real_21_Periphery.csv", package = "cplass")
#' traj <- read.csv(data_file)
#' path <- traj[traj$index_path == 3, ]
#'
#' fit <- CPLASS(path$t, path$x, path$y, iter_max = 2000, burn_in = 200,
#'                show_progress = FALSE)
#' plot_path_inferred(fit, motor = "Example", max_speed = 2, path_index = 3)
#' }
#' @export
plot_path_inferred <- function(path_info, xy_width = NA, t_lim = NA, motor, max_speed,
                                 state_shaded = TRUE, title_ind = TRUE,
                                 show_time_changes = FALSE, path_index = NULL) {
  if (is.na(xy_width) || is.na(t_lim[1])) {
    frame_info <- get_plot_size(list(path_info))
  }
  if (is.na(xy_width)) xy_width <- frame_info$xy_width
  if (is.na(t_lim[1])) t_lim <- c(frame_info$t_begin, frame_info$t_end)

  path <- path_info$path_inferred
  segments <- path_info$segments_inferred
  x_min <- mean(path$x) - xy_width / 2
  x_max <- mean(path$x) + xy_width / 2
  y_min <- mean(path$y) - xy_width / 2
  y_max <- mean(path$y) + xy_width / 2

  title_text <- if (title_ind) {
    paste0(motor, " Path ", .title_index(path_index))
  } else {
    paste0(motor, " Path ")
  }

  p_xy <- ggplot2::ggplot() +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$a, y = .data$b), col = col_inferred) +
    ggplot2::xlim(x_min, x_max) +
    ggplot2::ylim(y_min, y_max) +
    ggplot2::ggtitle(title_text) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  p_xt <- ggplot2::ggplot()
  p_yt <- ggplot2::ggplot()

  if (state_shaded) {
    shades <- c(col_stationary_inferred, col_motile_inferred)
    if (segments$cp_times[length(segments$cp_times)] > path$t[length(path$t)]) {
      segments$cp_times[length(segments$cp_times)] <- path$t[length(path$t)]
    }
    seg_ends <- c(path$t[1], segments$cp_times)
    for (i in seq_along(segments$cp_times)) {
      p_xt <- p_xt + ggplot2::annotate("rect",
        xmin = seg_ends[i], xmax = seg_ends[i + 1], ymin = x_min, ymax = x_max,
        alpha = 0.2, fill = shades[segments$states[i] + 1]
      )
      p_yt <- p_yt + ggplot2::annotate("rect",
        xmin = seg_ends[i], xmax = seg_ends[i + 1], ymin = y_min, ymax = y_max,
        alpha = 0.2, fill = shades[segments$states[i] + 1]
      )
    }
  }

  if (!show_time_changes) {
    p_xt <- p_xt +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$x)) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$a), col = col_inferred) +
      ggplot2::xlim(t_lim[1], t_lim[2]) +
      ggplot2::ylim(x_min, x_max) +
      ggplot2::ggtitle("x- and a- Time Series") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

    p_yt <- p_yt +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$y)) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$b), col = col_inferred) +
      ggplot2::xlim(t_lim[1], t_lim[2]) +
      ggplot2::ylim(y_min, y_max) +
      ggplot2::ggtitle("y- and b- Time Series") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else {
    cp_no_last <- segments$cp_times[-length(segments$cp_times)]
    p_xt <- p_xt +
      ggplot2::geom_vline(xintercept = cp_no_last, linetype = "twodash", col = "black", linewidth = 0.25) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$x)) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$a), col = col_inferred) +
      ggplot2::xlim(t_lim[1], t_lim[2]) +
      ggplot2::ylim(x_min, x_max) +
      ggplot2::ggtitle("x- and a- Time Series") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

    p_yt <- p_yt +
      ggplot2::geom_vline(xintercept = cp_no_last, linetype = "twodash", col = "black", linewidth = 0.25) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$y)) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$b), col = col_inferred) +
      ggplot2::xlim(t_lim[1], t_lim[2]) +
      ggplot2::ylim(y_min, y_max) +
      ggplot2::ggtitle("y- and b- Time Series") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  p_segments <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, col = "gray") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, col = "gray") +
    ggplot2::geom_point(data = segments, ggplot2::aes(x = .data$durations, y = .data$speeds),
                          alpha = 0.4, size = 0.6, col = col_inferred) +
    ggplot2::xlim(0, t_lim[2]) +
    ggplot2::ylim(0, max_speed) +
    ggplot2::xlab("Segment duration (s)") +
    ggplot2::ylab(expression(paste("Speed (", mu, "m/s)"))) +
    ggplot2::ggtitle(paste0("Segments: ", length(segments$cp_times))) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1, size = 10, margin = ggplot2::margin(0, 0, -10, 0)))

  layout <- c(
    patchwork::area(t = 1, l = 1, b = 5, r = 3),
    patchwork::area(t = 1, l = 4, b = 4, r = 7),
    patchwork::area(t = 5, l = 4, b = 7, r = 7),
    patchwork::area(t = 6, l = 1, b = 7, r = 3)
  )

  (p_xy + p_xt + p_yt + p_segments + patchwork::plot_layout(design = layout))
}

#' Alternate 3-row layout for a single segmented trajectory
#' @inheritParams plot_path_inferred
#' @export
plot_path_inferred_second_version <- function(path_info, xy_width = NA, t_lim = NA, motor, max_speed,
                                                state_shaded = TRUE, title_ind = TRUE,
                                                show_time_changes = FALSE, path_index = NULL) {
  if (is.na(xy_width) || is.na(t_lim[1])) {
    frame_info <- get_plot_size(list(path_info))
  }
  if (is.na(xy_width)) xy_width <- frame_info$xy_width
  if (is.na(t_lim[1])) t_lim <- c(frame_info$t_begin, frame_info$t_end)

  path <- path_info$path_inferred
  segments <- path_info$segments_inferred
  x_min <- min(path$x) - 0.01
  x_max <- max(path$x) + 0.01
  y_min <- min(path$y) - 0.01
  y_max <- max(path$y) + 0.01

  title_text <- if (title_ind) {
    paste0(motor, " Path ", .title_index(path_index))
  } else {
    paste0(motor, " Path ")
  }

  p_xy <- ggplot2::ggplot() +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$a, y = .data$b), col = col_inferred) +
    ggplot2::xlim(x_min, x_max) +
    ggplot2::ylim(y_min, y_max) +
    ggplot2::ggtitle(title_text) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  p_xt <- ggplot2::ggplot()
  p_yt <- ggplot2::ggplot()

  if (state_shaded) {
    shades <- c(col_stationary_inferred, col_motile_inferred)
    if (segments$cp_times[length(segments$cp_times)] > path$t[length(path$t)]) {
      segments$cp_times[length(segments$cp_times)] <- path$t[length(path$t)]
    }
    seg_ends <- c(path$t[1], segments$cp_times)
    for (i in seq_along(segments$cp_times)) {
      p_xt <- p_xt + ggplot2::annotate("rect",
        xmin = seg_ends[i], xmax = seg_ends[i + 1], ymin = x_min, ymax = x_max,
        alpha = 0.2, fill = shades[segments$states[i] + 1]
      )
      p_yt <- p_yt + ggplot2::annotate("rect",
        xmin = seg_ends[i], xmax = seg_ends[i + 1], ymin = y_min, ymax = y_max,
        alpha = 0.2, fill = shades[segments$states[i] + 1]
      )
    }
  }

  if (!show_time_changes) {
    p_xt <- p_xt +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$x)) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$a), col = col_inferred) +
      ggplot2::xlim(t_lim[1], t_lim[2]) +
      ggplot2::ylim(x_min, x_max) +
      ggplot2::ggtitle("x- and a- Time Series") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

    p_yt <- p_yt +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$y)) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$b), col = col_inferred) +
      ggplot2::xlim(t_lim[1], t_lim[2]) +
      ggplot2::ylim(y_min, y_max) +
      ggplot2::ggtitle("y- and b- Time Series") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else {
    cp_no_last <- segments$cp_times[-length(segments$cp_times)]
    p_xt <- p_xt +
      ggplot2::geom_vline(xintercept = cp_no_last, linetype = "twodash", col = "black", linewidth = 0.25) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$x)) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$a), col = col_inferred) +
      ggplot2::xlim(t_lim[1], t_lim[2]) +
      ggplot2::ylim(x_min, x_max) +
      ggplot2::ggtitle("x- and a- Time Series") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

    p_yt <- p_yt +
      ggplot2::geom_vline(xintercept = cp_no_last, linetype = "twodash", col = "black", linewidth = 0.25) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$y)) +
      ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$b), col = col_inferred) +
      ggplot2::xlim(t_lim[1], t_lim[2]) +
      ggplot2::ylim(y_min, y_max) +
      ggplot2::ggtitle("y- and b- Time Series") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  p_segments <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, col = "gray") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, col = "gray") +
    ggplot2::geom_point(data = segments, ggplot2::aes(x = .data$durations, y = .data$speeds),
                          alpha = 0.4, size = 0.6, col = col_inferred) +
    ggplot2::xlim(0, t_lim[2]) +
    ggplot2::ylim(0, max_speed) +
    ggplot2::xlab("Segment duration (s)") +
    ggplot2::ylab(expression(paste("Speed (", mu, "m/s)"))) +
    ggplot2::ggtitle(paste0("Segments: ", length(segments$cp_times))) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1, size = 10, margin = ggplot2::margin(0, 0, -10, 0)))

  layout <- c(
    patchwork::area(t = 1, l = 1, b = 1, r = 2),
    patchwork::area(t = 1, l = 3, b = 1, r = 4),
    patchwork::area(t = 2, l = 1, b = 2, r = 4),
    patchwork::area(t = 3, l = 1, b = 3, r = 4)
  )

  (p_xy + p_segments + p_xt + p_yt + patchwork::plot_layout(design = layout))
}

#' Just the x-y panel of a segmented trajectory
#' @inheritParams plot_path_inferred
#' @export
plot_path_inferred_xy <- function(path_info, xy_width = NA, t_lim = NA, motor, max_speed, state_shaded = TRUE) {
  if (is.na(xy_width) || is.na(t_lim[1])) {
    frame_info <- get_plot_size(list(path_info))
  }
  if (is.na(xy_width)) xy_width <- frame_info$xy_width
  if (is.na(t_lim[1])) t_lim <- c(frame_info$t_begin, frame_info$t_end)

  path <- path_info$path_inferred
  x_min <- mean(path$x) - xy_width / 2
  x_max <- mean(path$x) + xy_width / 2
  y_min <- mean(path$y) - xy_width / 2
  y_max <- mean(path$y) + xy_width / 2

  ggplot2::ggplot() +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$a, y = .data$b), col = col_inferred) +
    ggplot2::xlim(x_min, x_max) +
    ggplot2::ylim(y_min, y_max) +
    ggplot2::ggtitle(paste0(motor, " Path ")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Side-by-side comparison of true vs. inferred segmentation
#'
#' @param path_info a path object with both ground-truth (`path`,
#'   `segments`) and inferred (`path_inferred`, `segments_inferred`) fields.
#' @inheritParams plot_path_inferred
#' @export
plot_path_actual_and_inferred <- function(path_info, xy_width = NA, t_lim = NA,
                                            state_shaded = TRUE, path_index = NULL) {
  if (is.na(xy_width) || is.na(t_lim[1])) {
    frame_info <- get_plot_size(path_info)
  }
  if (is.na(xy_width)) xy_width <- frame_info$xy_width
  if (is.na(t_lim[1])) t_lim <- c(frame_info$t_begin, frame_info$t_end)

  path <- path_info$path
  segments <- path_info$segments
  path_inf <- path_info$path_inferred
  segments_inf <- path_info$segments_inferred

  duration_difference <- duration_of_difference(path_info)

  x_min <- mean(path$x) - xy_width / 2
  x_max <- mean(path$x) + xy_width / 2
  y_min <- mean(path$y) - xy_width / 2
  y_max <- mean(path$y) + xy_width / 2

  p_xy <- ggplot2::ggplot() +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$a, y = .data$b), col = col_actual) +
    ggplot2::geom_path(data = path_inf, ggplot2::aes(x = .data$a, y = .data$b), col = col_inferred) +
    ggplot2::xlim(x_min, x_max) +
    ggplot2::ylim(y_min, y_max) +
    ggplot2::ggtitle(paste0("Simulated path ", .title_index(path_index))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  p_xt <- ggplot2::ggplot()
  p_yt <- ggplot2::ggplot()

  if (state_shaded) {
    shades <- c(col_stationary, col_motile)
    if (segments$cp_times[length(segments$cp_times)] > path$t[length(path$t)]) {
      segments$cp_times[length(segments$cp_times)] <- path$t[length(path$t)]
    }
    seg_ends <- c(path$t[1], segments$cp_times)
    for (i in seq_along(segments$cp_times)) {
      p_xt <- p_xt + ggplot2::annotate("rect",
        xmin = seg_ends[i], xmax = seg_ends[i + 1], ymin = x_min, ymax = x_max,
        alpha = 0.2, fill = shades[segments$states[i] + 1]
      )
      p_yt <- p_yt + ggplot2::annotate("rect",
        xmin = seg_ends[i], xmax = seg_ends[i + 1], ymin = y_min, ymax = y_max,
        alpha = 0.2, fill = shades[segments$states[i] + 1]
      )
    }

    shades <- c(col_stationary_inferred, col_motile_inferred)
    if (segments$cp_times[length(segments$cp_times)] > path$t[length(path$t)]) {
      segments$cp_times[length(segments$cp_times)] <- path$t[length(path$t)]
    }
    seg_ends <- c(path_inf$t[1], segments_inf$cp_times)
    for (i in seq_along(segments_inf$cp_times)) {
      p_xt <- p_xt + ggplot2::annotate("rect",
        xmin = seg_ends[i], xmax = seg_ends[i + 1], ymin = x_min, ymax = x_max,
        alpha = 0.2, fill = shades[segments_inf$states[i] + 1]
      )
      p_yt <- p_yt + ggplot2::annotate("rect",
        xmin = seg_ends[i], xmax = seg_ends[i + 1], ymin = y_min, ymax = y_max,
        alpha = 0.2, fill = shades[segments_inf$states[i] + 1]
      )
    }
  }

  p_xt <- p_xt +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$x)) +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$a), col = col_actual) +
    ggplot2::geom_path(data = path_inf, ggplot2::aes(x = .data$t, y = .data$a), col = col_inferred) +
    ggplot2::xlim(t_lim[1], t_lim[2]) +
    ggplot2::ylim(x_min, x_max) +
    ggplot2::ggtitle("x- and a- Time Series") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  p_yt <- p_yt +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$y)) +
    ggplot2::geom_path(data = path, ggplot2::aes(x = .data$t, y = .data$b), col = col_actual) +
    ggplot2::geom_path(data = path_inf, ggplot2::aes(x = .data$t, y = .data$b), col = col_inferred) +
    ggplot2::xlim(t_lim[1], t_lim[2]) +
    ggplot2::ylim(y_min, y_max) +
    ggplot2::ggtitle("y- and b- Time Series") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  p_segments <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, col = "gray") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, col = "gray") +
    ggplot2::geom_point(data = segments, ggplot2::aes(x = .data$durations, y = .data$speeds),
                          alpha = 0.4, size = 0.6, col = col_actual) +
    ggplot2::geom_point(data = segments_inf, ggplot2::aes(x = .data$durations, y = .data$speeds),
                          alpha = 0.4, size = 0.6, col = col_inferred) +
    ggplot2::xlim(t_lim[1], t_lim[2]) +
    ggplot2::ylim(0, 1) +
    ggplot2::xlab("Segment duration (s)") +
    ggplot2::ylab(expression(paste("Speed (", mu, "m/s)"))) +
    ggplot2::ggtitle(paste0(
      "Segments: ", length(segments$cp_times),
      "\n Inf. Segments: ", length(segments_inf$cp_times),
      "\n inference Gap: ", round(duration_difference, digits = 2)
    )) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1, size = 10, margin = ggplot2::margin(0, 0, -10, 0)))

  layout <- c(
    patchwork::area(t = 1, l = 1, b = 5, r = 3),
    patchwork::area(t = 1, l = 4, b = 4, r = 7),
    patchwork::area(t = 5, l = 4, b = 7, r = 7),
    patchwork::area(t = 6, l = 1, b = 7, r = 3)
  )

  (p_xy + p_xt + p_yt + p_segments + patchwork::plot_layout(design = layout))
}

#' Plot Cumulative Speed Allocation curves
#'
#' @param csa a tibble with columns `s`, `csa`, `label`, `subsample`, `Hz`.
#' @param max_error unused, kept for interface compatibility.
#' @param legend logical; overlay the "True" and first subsample curves
#'   with a legend. Default `TRUE`.
#' @export
plot_csa <- function(csa, max_error = 0.1, legend = TRUE) {
  pl_title <- paste("Cumulative Speed Allocation:", csa$Hz[which(!is.na(csa$Hz))][1], "Hz")
  p_csa <- ggplot2::ggplot()

  if (legend) {
    csa_three <- dplyr::filter(csa, .data$label == "True" & .data$s > 0)
    csa_three <- dplyr::bind_rows(csa_three, dplyr::filter(csa, .data$s > 0 & .data$subsample == 1))
    p_csa <- p_csa +
      ggplot2::geom_path(data = csa_three, ggplot2::aes(x = .data$s, y = .data$csa, group = .data$label, col = .data$label))
  }

  p_csa +
    ggplot2::geom_path(
      data = dplyr::filter(csa, .data$label == "Sample" & .data$s > 0),
      ggplot2::aes(x = .data$s, y = .data$csa, group = .data$subsample), col = col_actual, alpha = 0.25
    ) +
    ggplot2::geom_path(
      data = dplyr::filter(csa, .data$label == "Inferred" & .data$s > 0),
      ggplot2::aes(x = .data$s, y = .data$csa, group = .data$subsample), col = col_inferred, alpha = 0.25
    ) +
    ggplot2::geom_path(
      data = dplyr::filter(csa, .data$label == "True" & .data$s > 0),
      ggplot2::aes(x = .data$s, y = .data$csa), col = "black"
    ) +
    ggplot2::ggtitle(pl_title) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = c(0.9, 0.2))
}

#' Total duration over which true and inferred segmentations disagree on state
#'
#' Sweeps two segmentations (true and CPLASS-inferred) simultaneously with a
#' two-pointer merge over their changepoints and sums the overlap duration
#' where the two disagree on motile/stationary state.
#'
#' @param path_info a path object with `path`, `segments`,
#'   `segments_inferred`.
#' @return numeric, total duration of disagreement (same units as `t`).
#' @export
duration_of_difference <- function(path_info) {
  cp_times <- path_info$segments$cp_times
  cp_times[length(cp_times)] <- max(path_info$path$t)
  segments1 <- data.frame(
    start = c(path_info$path$t[1], cp_times[-length(cp_times)]),
    end = cp_times,
    state = path_info$segments$states
  )
  cp_times <- path_info$segments_inferred$cp_times
  segments2 <- data.frame(
    start = c(path_info$path$t[1], cp_times[-length(cp_times)]),
    end = cp_times,
    state = path_info$segments_inferred$states
  )

  total_overlap_length <- 0
  i <- 1
  j <- 1

  while (i <= nrow(segments1) && j <= nrow(segments2)) {
    count1 <- 1
    count2 <- 1
    chain1 <- segments1[i, ]
    chain2 <- segments2[j, ]

    if (chain1$end == chain2$end) {
      if (chain1$state != chain2$state) {
        overlap_length <- min(chain1$end, chain2$end) - max(chain1$start, chain2$start)
        total_overlap_length <- total_overlap_length + overlap_length
      }
      i <- i + 1
      j <- j + 1
    } else if (chain1$end < chain2$end) {
      if (chain1$state != chain2$state) {
        overlap_length <- min(chain1$end, chain2$end) - max(chain1$start, chain2$start)
        total_overlap_length <- total_overlap_length + overlap_length
      }
      next_chain1 <- segments1[i + count1, ]
      while (next_chain1$end <= chain2$end) {
        if (next_chain1$state != chain2$state) {
          overlap_length <- next_chain1$end - next_chain1$start
          total_overlap_length <- total_overlap_length + overlap_length
        }
        count1 <- count1 + 1
        if ((i + count1) > nrow(segments1)) break
        next_chain1 <- segments1[i + count1, ]
      }
      if ((i + count1) > nrow(segments1)) {
        break
      } else if (next_chain1$start <= chain2$end) {
        if (next_chain1$state != chain2$state) {
          overlap_length <- chain2$end - next_chain1$start
          total_overlap_length <- total_overlap_length + overlap_length
        }
      }
      i <- i + count1
      j <- j + 1
    } else {
      # chain1$end > chain2$end
      if (chain1$state != chain2$state) {
        overlap_length <- min(chain1$end, chain2$end) - max(chain1$start, chain2$start)
        total_overlap_length <- total_overlap_length + overlap_length
      }
      next_chain2 <- segments2[j + count2, ]
      while (next_chain2$end <= chain1$end) {
        if (next_chain2$state != chain1$state) {
          overlap_length <- next_chain2$end - next_chain2$start
          total_overlap_length <- total_overlap_length + overlap_length
        }
        count2 <- count2 + 1
        if ((j + count2) > nrow(segments2)) break
        next_chain2 <- segments2[j + count2, ]
      }
      if ((j + count2) > nrow(segments2)) {
        break
      } else if (next_chain2$start <= chain1$end) {
        if (next_chain2$state != chain1$state) {
          overlap_length <- chain1$end - next_chain2$start
          total_overlap_length <- total_overlap_length + overlap_length
        }
      }
      j <- j + count2
      i <- i + 1
    }
  }

  total_overlap_length
}
