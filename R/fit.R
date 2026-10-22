#' Fit a continuous piecewise-linear anchor path given changepoints
#'
#' Internal, performance-critical fitting routine. Given a fixed set of
#' changepoint indices, this constructs the hinge-basis design matrix from
#' Eq. (2.10) of the companion paper and computes the maximum-likelihood
#' piecewise-linear fit in both coordinates from a *single* QR
#' decomposition (Eq. 2.13-2.14), rather than one decomposition per
#' coordinate.
#'
#' This function does the same linear algebra as the original
#' \code{piecewise_linear_con()} but:
#' \itemize{
#'   \item builds the hinge design matrix with a vectorized \code{outer()}
#'     call instead of a nested \code{for} loop over \code{i, j},
#'   \item calls \code{qr()} once and reuses it for both the x- and
#'     y-coordinate right-hand sides via \code{qr.coef(qrA, cbind(x, y))},
#'     instead of two independent \code{qr.solve()} calls,
#'   \item returns plain vectors/matrices rather than a `tibble`, since this
#'     is called many times per MCMC iteration.
#' }
#'
#' @param t numeric vector of observation times.
#' @param x numeric vector of x-coordinates.
#' @param y numeric vector of y-coordinates.
#' @param cp integer vector of changepoint indices (may be empty/length 0).
#' @return A list with entries `logical` (fit succeeded), `alpha`, `beta`
#'   (regression coefficients for x and y), `vx`, `vy` (segment velocity
#'   components), `seg_speed`, `seg_time`, `seg_theta`, `RSS`, `sigma_hat`,
#'   `x_fit`, `y_fit`.
#' @keywords internal
#' @noRd
.fit_pla <- function(t, x, y, cp) {
  n <- length(t)
  cp <- sort(cp)
  m <- length(cp) + 1L # number of segments

  if (length(cp) > 0L) {
    # hinge basis columns: (t_i - t_cp)_+ , vectorized over all changepoints
    # at once. Equivalent to the original (t[i]-t[cp])*1{i>cp}: since t is
    # strictly increasing on the observation grid, t[i]-t[cp] > 0 iff i > cp,
    # so pmax(., 0) reproduces the same values exactly.
    hinge <- outer(t, t[cp], function(ti, tc) pmax(ti - tc, 0))
    A <- cbind(1, t, hinge)
  } else {
    A <- cbind(1, t)
  }

  qrA <- tryCatch(qr(A), error = function(e) NULL)
  if (is.null(qrA) || qrA$rank < ncol(A)) {
    return(list(logical = FALSE))
  }

  coefs <- tryCatch(qr.coef(qrA, cbind(x, y)), error = function(e) NULL)
  if (is.null(coefs) || anyNA(coefs)) {
    return(list(logical = FALSE))
  }

  alpha <- coefs[, 1]
  beta <- coefs[, 2]

  fitted <- A %*% coefs
  x_fit <- fitted[, 1]
  y_fit <- fitted[, 2]

  resid <- cbind(x, y) - fitted
  RSS <- sum(resid * resid)

  # Reject any fit whose residual degrees of freedom are not strictly
  # positive, or whose RSS is non-finite/non-positive. With m segments fit
  # jointly across both coordinates, the residual degrees of freedom are
  # 2*n - 2*m - 2 = 2*(n - c(r) - 2), which is strictly positive exactly
  # when c(r) <= n - 3. The existing rank check above catches n < ncol(A)
  # (rank-deficient); it does NOT catch the case n == ncol(A) exactly,
  # where the design is square and full rank but interpolates the data
  # exactly (RSS = 0, resid_df = 0), which is an unbounded-criterion
  # failure mode. This check is an independent safety net that also
  # protects any caller (e.g. a direct call to piecewise_linear_con())
  # that bypasses the upstream .r_is_admissible() guard in CS(). RSS is
  # compared to 0 rather than a numerical tolerance, per the mathematical
  # requirement RSS > 0; the is.finite() checks separately guard against
  # NaN/Inf from ill-conditioned but nominally full-rank fits.
  resid_df <- 2 * n - 2 * m - 2
  if (!is.finite(RSS) || RSS <= 0 || !is.finite(resid_df) || resid_df <= 0) {
    return(list(logical = FALSE))
  }

  if (length(cp) > 0L) {
    sigma_hat <- sqrt(RSS / (2 * n - 2 * m - 2))
    vx <- cumsum(alpha[-1]) # alpha[2:(m+1)], cumulative sums per segment
    vy <- cumsum(beta[-1])
    seg_speed <- sqrt(vx * vx + vy * vy)
    seg_time <- diff(c(t[1], t[cp], t[n]))
    seg_theta <- atan2(vy, vx)
  } else {
    sigma_hat <- sqrt(RSS / (2 * n - 4))
    vx <- alpha[2]
    vy <- beta[2]
    seg_speed <- sqrt(vx * vx + vy * vy)
    seg_time <- t[n] - t[1]
    seg_theta <- atan2(vy, vx)
  }

  list(
    logical = TRUE,
    alpha = alpha,
    beta = beta,
    vx = vx,
    vy = vy,
    seg_speed = seg_speed,
    seg_time = seg_time,
    seg_theta = seg_theta,
    RSS = RSS,
    sigma_hat = sigma_hat,
    x_fit = x_fit,
    y_fit = y_fit
  )
}

#' Continuous piecewise-linear fit given a set of changepoints
#'
#' Given observed times \code{t} and coordinates \code{x}, \code{y}, and a
#' vector of changepoint indices \code{this_cp}, compute the maximum
#' likelihood continuous piecewise-linear anchor path (Section 2.2 of the
#' companion paper).
#'
#' This is a drop-in replacement for the original research-script function
#' of the same name: the argument list and the structure of the returned
#' list are unchanged, so existing scripts that call
#' \code{piecewise_linear_con()} do not need to be modified. Internally it
#' now delegates to an internal fast fitting routine, which fits both coordinates from
#' a single QR decomposition instead of two.
#'
#' @param t numeric vector of observation times.
#' @param x numeric vector of x-coordinates.
#' @param y numeric vector of y-coordinates.
#' @param this_cp integer vector of changepoint indices.
#' @param old_version logical. If \code{TRUE} (default), return the flat
#'   list format used internally by \code{\link{CS}} and \code{\link{MHsearch}}.
#'   If \code{FALSE}, return the `tibble`-based format
#'   (`path_inferred`, `segments_inferred`) used by the plotting and CSA
#'   helpers.
#' @return See Details; matches the original function's return format.
#' @export
piecewise_linear_con <- function(t, x, y, this_cp, old_version = TRUE) {
  n <- length(x)
  fit <- .fit_pla(t, x, y, this_cp)

  if (!isTRUE(fit$logical)) {
    return(list(logical = FALSE))
  }

  this_cp <- sort(this_cp)

  if (old_version) {
    list(
      t = t,
      x = x,
      y = y,
      x_piecewise = fit$x_fit,
      y_piecewise = fit$y_fit,
      this_cp = this_cp,
      path_segvel = fit$seg_speed,
      path_segtime = fit$seg_time,
      path_segtheta = fit$seg_theta,
      path_segeta = 1 / fit$sigma_hat, # eta := 1/sigma_hat, matches original
      inter_alpha = fit$alpha[1],
      inter_beta = fit$beta[1],
      u_x = fit$vx,
      v_y = fit$vy,
      path_RSS = fit$RSS,
      logical = TRUE
    )
  } else {
    path_inferred <- tibble::tibble(
      t = t,
      j = rep("NA", n),
      x = x,
      y = y,
      a = fit$x_fit,
      b = fit$y_fit
    )

    segments_inferred <- tibble::tibble(
      cp_times = c(t[this_cp], t[n]),
      durations = fit$seg_time,
      states = infer_states_speed_cutoff(fit$seg_speed, cutoff = 0.1),
      speeds = fit$seg_speed,
      vx = fit$vx,
      vy = fit$vy,
      angles = fit$seg_theta
    )

    path_inferred$j <- build_j_inferred(path_inferred, segments_inferred)

    list(
      segments_inferred = segments_inferred,
      path_inferred = path_inferred
    )
  }
}
