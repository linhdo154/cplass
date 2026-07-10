#' Fill gaps in an unevenly-sampled 2D trajectory
#'
#' Regularizes a trajectory onto a uniform time grid, linearly interpolating
#' missing observations (with added noise matched to the observed
#' increment scale) and replacing implausible jumps.
#'
#' @param t,x,y observed times and coordinates (may have gaps/dropped frames).
#' @param track_id an identifier carried through to the output.
#' @return a list with the regularized `t`, `x`, `y`, indices that were
#'   `filled` vs. `replaced` for each coordinate, and `track_id`.
#' @export
fill_missing_data <- function(t, x, y, track_id) {
  t <- round(t, digits = 4)
  this_dt <- round(min(diff(t)), digits = 4)
  this_new_t <- seq(t[1], t[length(t)], by = this_dt)
  this_new_x <- numeric(length(this_new_t))
  this_new_y <- numeric(length(this_new_t))

  sigma_est_x <- stats::sd(diff(x)) / 2
  sigma_est_y <- stats::sd(diff(y)) / 2

  find_next_good <- function(v, idx) {
    i <- idx
    while (v[i] == 0) i <- i + 1
    i
  }

  for (tt in seq_along(t)) {
    this_new_idx <- which(abs(this_new_t - t[tt]) < 0.0005)
    this_new_x[this_new_idx] <- x[tt]
    this_new_y[this_new_idx] <- y[tt]
  }

  fill_idx_x <- which(this_new_x == 0)
  fill_idx_y <- which(this_new_y == 0)

  last_good_idx <- 1
  for (tt in seq(2, length(this_new_t) - 1)) {
    if (this_new_x[tt] == 0 || this_new_y[tt] == 0) {
      next_good_idx_x <- find_next_good(this_new_x, tt)
      next_good_idx_y <- find_next_good(this_new_y, tt)

      this_new_x[tt] <- (this_new_t[next_good_idx_x] - this_new_t[tt]) /
        (this_new_t[next_good_idx_x] - this_new_t[last_good_idx]) * this_new_x[last_good_idx] +
        (this_new_t[tt] - this_new_t[last_good_idx]) /
        (this_new_t[next_good_idx_x] - this_new_t[last_good_idx]) * this_new_x[next_good_idx_x] +
        sigma_est_x * stats::rnorm(1, 0, 1)

      this_new_y[tt] <- (this_new_t[next_good_idx_y] - this_new_t[tt]) /
        (this_new_t[next_good_idx_y] - this_new_t[last_good_idx]) * this_new_y[last_good_idx] +
        (this_new_t[tt] - this_new_t[last_good_idx]) /
        (this_new_t[next_good_idx_y] - this_new_t[last_good_idx]) * this_new_y[next_good_idx_y] +
        sigma_est_y * stats::rnorm(1, 0, 1)
    } else {
      last_good_idx <- tt
    }
  }

  replace_idx_x <- which(diff(this_new_x) > 5 * sigma_est_x) + 1
  replace_idx_y <- which(diff(this_new_y) > 5 * sigma_est_y) + 1

  for (tt in replace_idx_x) {
    if ((tt > 1 && tt < length(this_new_t)) || tt < length(this_new_x)) {
      last_good_idx_x <- tt - 1
      next_good_idx_x <- tt + 1
      this_new_x[tt] <- (this_new_t[next_good_idx_x] - this_new_t[tt]) /
        (this_new_t[next_good_idx_x] - this_new_t[last_good_idx_x]) * this_new_x[last_good_idx_x] +
        (this_new_t[tt] - this_new_t[last_good_idx_x]) /
        (this_new_t[next_good_idx_x] - this_new_t[last_good_idx_x]) * this_new_x[next_good_idx_x] +
        sigma_est_x * stats::rnorm(1, 0, 1)
    }
  }

  for (tt in replace_idx_y) {
    if ((tt > 1 && tt < length(this_new_t)) || tt < length(this_new_y)) {
      last_good_idx_y <- tt - 1
      next_good_idx_y <- tt + 1
      this_new_y[tt] <- (this_new_t[next_good_idx_y] - this_new_t[tt]) /
        (this_new_t[next_good_idx_y] - this_new_t[last_good_idx_y]) * this_new_y[last_good_idx_y] +
        (this_new_t[tt] - this_new_t[last_good_idx_y]) /
        (this_new_t[next_good_idx_y] - this_new_t[last_good_idx_y]) * this_new_y[next_good_idx_y] +
        sigma_est_y * stats::rnorm(1, 0, 1)
    }
  }

  list(
    t = this_new_t, x = this_new_x, y = this_new_y,
    filled_x = fill_idx_x, filled_y = fill_idx_y,
    replaced_x = replace_idx_x, replaced_y = replace_idx_y,
    track_id = track_id
  )
}
