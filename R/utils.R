#' Classify inferred segments as motile/stationary by a speed cutoff
#'
#' @param seg_speeds numeric vector of inferred segment speeds.
#' @param cutoff speed threshold; segments faster than this are labeled 1
#'   (motile), others 0 (stationary). Default 0.1.
#' @return integer vector of 0/1 labels, same length as `seg_speeds`.
#' @export
infer_states_speed_cutoff <- function(seg_speeds, cutoff = 0.1) {
  as.integer(seg_speeds > cutoff)
}

#' Broadcast per-segment states back onto the original observation grid
#'
#' @param path_inf a `path_inferred` tibble (must have column `t`).
#' @param seg_inf a `segments_inferred` tibble (must have columns
#'   `cp_times`, `states`).
#' @return integer vector of per-observation state labels, length
#'   `nrow(path_inf)`.
#' @export
build_j_inferred <- function(path_inf, seg_inf) {
  j <- integer(nrow(path_inf))
  cpt <- c(path_inf$t[1] - 0.001, seg_inf$cp_times)
  for (i in seq_along(cpt)[-length(cpt)]) {
    seg_idx <- which(path_inf$t > cpt[i] & path_inf$t <= cpt[i + 1])
    j[seg_idx] <- seg_inf$states[i]
  }
  j
}
