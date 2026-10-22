#' Gaussian log-likelihood from residual sum of squares
#'
#' @param n number of observations.
#' @param RSS residual sum of squares (summed over both coordinates).
#' @param sd optional known noise standard deviation. If \code{NA} (default),
#'   the noise variance is profiled out (concentrated log-likelihood), which
#'   is the branch used everywhere in \code{\link{CS}}.
#' @return numeric log-likelihood value.
#' @export
loglikelihood <- function(n, RSS, sd = NA) {
  if (!is.na(sd)) {
    # Known-sd branch. NOTE: this branch is not exercised anywhere in the
    # package (CS() always calls with sd = NA), and the original script had
    # a typo here (`log*(2*pi)`, which is not valid R and would error if
    # ever reached). This is a best-effort fix to a standard bivariate
    # Gaussian log-likelihood with known sd; please confirm this matches
    # intent before relying on it.
    l <- -n * log(2 * pi) - n * log(sd^2) - RSS / (2 * sd^2)
  } else {
    l <- -n * log(2 * pi) - n * log(RSS) + n * log(2 * n) - n
  }
  l
}

#' Corrected AIC penalty
#' @param ncp number of changepoints.
#' @param n number of observations.
#' @examples
#' AICc(ncp = 2, n = 200)
#' @export
AICc <- function(ncp, n) {
  k <- 3*ncp + 5
  denom <- 2 * n - k - 1
  if (denom < 0) {
    return(0)
  }
  2 * k + 2 * k * (k + 1) / denom
}

#' Strengthened Schwarz Information Criterion penalty
#'
#' \eqn{(\log n)^\gamma \rho}, as in Definition 2 of the companion paper.
#' @param ncp number of changepoints.
#' @param n number of observations.
#' @param gamma exponent, \code{gamma > 1} recommended (default in
#'   \code{\link{CPLASS}} is 1.01).
#' @examples
#' sSIC(ncp = 2, n = 200, gamma = 1.01)
#' @export
sSIC <- function(ncp, n, gamma) {
  k <- 3*ncp + 5
  log(n)^gamma * k
}

#' Criterion function for a changepoint configuration
#'
#' Computes \eqn{\Phi(r) = 2 \hat L_n - \mathrm{pen}(r)} (Definition 1/2 of
#' the companion paper): the penalized log-likelihood of the continuous
#' piecewise-linear fit associated with a changepoint vector \code{r}, using
#' either the strengthened SIC penalty, AICc, or the max of the two
#' ("hybrid"), plus an optional dimensionless speed penalty.
#'
#' @param t,x,y trajectory data.
#' @param r changepoint indicator vector (0/1), length \code{length(t) - 2}.
#' @param p1 penalty coefficient \eqn{(\log n)^\gamma} for the linear
#'   (sSIC) term; only used when \code{pen == "ssic"} indirectly via the
#'   caller's precomputed value (kept for interface compatibility with the
#'   original script; \code{gamma} is what actually parametrizes sSIC here).
#' @param s_cap maximum segment speed with no penalty.
#' @param with_speed 1 to activate the speed penalty, 0 to deactivate.
#' @param eta non-negative, dimensionless strength of the speed penalty.
#'   The speed contribution is
#'   \eqn{\eta \sum_j (\hat s_j/s_{cap} - 1)_+}; \code{eta = 0}
#'   recovers the criterion without speed regularization.
#' @param sd optional known noise sd, passed to \code{\link{loglikelihood}}.
#' @param pen one of \code{"ssic"}, \code{"aicc"}, \code{"hybrid"}.
#' @param gamma exponent for the sSIC penalty.
#' @param K_max optional maximum number of segments \eqn{K(r) = c(r)+1}
#'   (the manuscript's practitioner-specified \eqn{\bar k} in Eq. 2.6). If
#'   \code{NULL} (default), no additional practitioner bound is enforced
#'   beyond the structural bound \eqn{K(r) \le n-2} implied by
#'   \eqn{c(r) \le n-3}. See \code{.r_is_admissible()} -- the manuscript
#'   intentionally leaves \eqn{\bar k} to the practitioner.
#' @param min_gap minimum allowed spacing between adjacent segment
#'   boundaries, in observation-index units. Default \code{1L} is the
#'   weakest possible requirement (structurally non-binding).
#' @return A list with `s` (criterion value), `llh`, `p` (penalty), `pv`
#'   (speed-penalty component), `logical`.
#' @export
CS <- function(t, x, y, r, p1, s_cap, with_speed, eta = 1, sd = NA,
               pen = "ssic", gamma, K_max = NULL, min_gap = 1L) {
  n <- length(t)

  if (!is.numeric(eta) || length(eta) != 1L || !is.finite(eta) || eta < 0) {
    stop("eta must be one finite non-negative number.", call. = FALSE)
  }
  if (isTRUE(with_speed != 0) &&
      (!is.numeric(s_cap) || length(s_cap) != 1L ||
       !is.finite(s_cap) || s_cap <= 0)) {
    stop("s_cap must be one finite positive number when the speed penalty is active.",
         call. = FALSE)
  }

  if (!.r_is_admissible(r, n, K_max = K_max, min_gap = min_gap)) {
    return(list(logical = FALSE))
  }

  cps <- which(r == 1) + 1
  pla <- piecewise_linear_con(t, x, y, cps)

  if (!isTRUE(pla$logical)) {
    return(list(logical = FALSE))
  }

  llh <- loglikelihood(n, pla$path_RSS, sd = sd)
  seg_vel <- pla$path_segvel
  pv <- if (isTRUE(with_speed != 0)) {
    eta * sum(pmax(0, seg_vel / s_cap - 1))
  } else {
    0
  }
  ncp <- length(cps)

  p <- switch(pen,
    hybrid = max(AICc(ncp, n), sSIC(ncp, n, gamma)) + with_speed * pv,
    aicc   = AICc(ncp, n) + with_speed * pv,
    ssic   = sSIC(ncp, n, gamma) + with_speed * pv,
    stop("Unknown `pen`: ", pen, call. = FALSE)
  )

  l <- 2*llh - p

  list(s = l, llh = llh, p = p, pv = pv, logical = TRUE, pla = pla)
}

#' Convert a changepoint index vector to a 0/1 indicator vector
#' @param cp integer vector of changepoint indices.
#' @param n number of observations.
#' @examples
#' cp_to_r(cp = c(5, 12), n = 20)
#' @export
cp_to_r <- function(cp, n) {
  output <- integer(n - 2)
  output[cp - 1] <- 1L
  output
}

#' Convert a 0/1 changepoint indicator vector to changepoint indices
#' @param r indicator vector.
#' @examples
#' r_to_cp(cp_to_r(cp = c(5, 12), n = 20))
#' @export
r_to_cp <- function(r) {
  1 + which(r == 1)
}

#' Check whether a changepoint indicator vector is structurally admissible
#'
#' Central admissibility check used throughout the package (by
#' \code{\link{CS}} and by the proposal mechanism, \code{\link{proposals}}).
#' A changepoint vector \eqn{r \in \{0,1\}^{n-2}} is structurally admissible
#' if all of the following hold, where \eqn{c(r) = \sum_i r_i} is the
#' number of changepoints and \eqn{K(r) = c(r)+1} is the number of
#' segments:
#' \enumerate{
#'   \item \code{r} has the correct length, \eqn{n-2};
#'   \item every entry of \code{r} is 0 or 1;
#'   \item \eqn{c(r) \le n-3}, which is exactly the condition for the
#'     residual degrees of freedom \eqn{2(n-c(r)-2)} of the continuous
#'     piecewise-linear fit (see \code{.fit_pla()}) to be strictly
#'     positive -- i.e. the hinge design is not saturated;
#'   \item \eqn{K(r)} does not exceed the effective ceiling
#'     \eqn{\min(K_{\max}, n-2)} (or just \eqn{n-2} if \code{K_max} is not
#'     supplied -- see "Note" below);
#'   \item every gap between adjacent segment boundaries (in observation-
#'     index units, including the trajectory endpoints) is at least
#'     \code{min_gap}.
#' }
#' This function performs only these structural/combinatorial checks; it
#' does not fit the piecewise-linear model and cannot detect numerical
#' issues (e.g. near-collinear but technically full-rank designs). The
#' independent numerical safety net in \code{.fit_pla()} (full-rank check,
#' plus residual-degrees-of-freedom and RSS checks) still applies and
#' should not be removed or relied upon less because of this helper.
#'
#' @note \code{K_max} corresponds to the manuscript's practitioner-
#'   specified upper bound \eqn{\bar k} on the number of segments (Eq.
#'   2.6), which the manuscript deliberately leaves unspecified as a fixed
#'   numerical value -- it is the practitioner's choice for a given
#'   analysis. \code{K_max = NULL} (default) means no additional
#'   practitioner bound is imposed beyond the structural ceiling
#'   \eqn{K(r) \le n-2}, which is already implied by check 3 above.
#'
#' @param r changepoint indicator vector (0/1), intended length
#'   \code{n - 2}.
#' @param n number of observations in the trajectory (i.e. \code{length(t)}
#'   for the corresponding \code{t}).
#' @param K_max optional maximum number of segments; see "Note".
#' @param min_gap minimum allowed spacing between adjacent segment
#'   boundaries, in observation-index units. Default \code{1L} (weakest
#'   possible requirement, structurally non-binding).
#' @return \code{TRUE} if \code{r} is structurally admissible, \code{FALSE}
#'   otherwise. Never errors on malformed input (e.g. wrong length,
#'   non-binary entries) -- returns \code{FALSE} instead, since this is
#'   meant to be usable as a cheap guard inside hot loops.
#' @keywords internal
#' @noRd
.r_is_admissible <- function(r, n, K_max = NULL, min_gap = 1L) {
  if (!is.numeric(r) || length(r) != n - 2L) {
    return(FALSE)
  }
  if (anyNA(r) || !all(r %in% c(0, 1))) {
    return(FALSE)
  }

  c_r <- sum(r)
  if (c_r > n - 3L) {
    return(FALSE)
  }

  K_r <- c_r + 1L
  K_max_effective <- if (is.null(K_max)) n - 2L else min(K_max, n - 2L)
  if (K_r > K_max_effective) {
    return(FALSE)
  }

  cps <- r_to_cp(r)
  boundaries <- c(1L, cps, n)
  if (any(diff(boundaries) < min_gap)) {
    return(FALSE)
  }

  TRUE
}

#' Validate the \code{K_max}/\code{min_gap} controls
#'
#' Lightweight, one-time validation for the public \code{K_max} and
#' \code{min_gap} controls, intended to be called once at a public entry
#' point (\code{\link{MHsearch}}) rather than repeatedly inside hot-loop
#' helpers such as \code{.r_is_admissible()}. Errors clearly on malformed
#' input rather than letting it silently propagate into confusing
#' downstream behavior.
#'
#' @param K_max the value to validate: \code{NULL}, or a single finite
#'   integer-like value \code{>= 1}.
#' @param min_gap the value to validate: a single finite integer-like
#'   value \code{>= 1}.
#' @return Invisibly \code{TRUE} if both are valid; otherwise throws an
#'   error via \code{stop()}.
#' @keywords internal
#' @noRd
.validate_K_max_min_gap <- function(K_max, min_gap) {
  if (!is.null(K_max)) {
    if (!is.numeric(K_max) || length(K_max) != 1L || is.na(K_max) ||
        !is.finite(K_max) || K_max != round(K_max) || K_max < 1) {
      stop(
        "`K_max` must be NULL or a single finite integer-like value >= 1 ",
        "(got: ", paste(deparse(K_max), collapse = " "), ").",
        call. = FALSE
      )
    }
  }
  if (!is.numeric(min_gap) || length(min_gap) != 1L || is.na(min_gap) ||
      !is.finite(min_gap) || min_gap != round(min_gap) || min_gap < 1) {
    stop(
      "`min_gap` must be a single finite integer-like value >= 1 ",
      "(got: ", paste(deparse(min_gap), collapse = " "), ").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
