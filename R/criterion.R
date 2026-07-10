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
#' ("hybrid"), plus an optional speed penalty.
#'
#' @param t,x,y trajectory data.
#' @param r changepoint indicator vector (0/1), length \code{length(t) - 2}.
#' @param p1 penalty coefficient \eqn{(\log n)^\gamma} for the linear
#'   (sSIC) term; only used when \code{pen == "ssic"} indirectly via the
#'   caller's precomputed value (kept for interface compatibility with the
#'   original script; \code{gamma} is what actually parametrizes sSIC here).
#' @param s_cap maximum segment speed with no penalty.
#' @param with_speed 1 to activate the speed penalty, 0 to deactivate.
#' @param sd optional known noise sd, passed to \code{\link{loglikelihood}}.
#' @param pen one of \code{"ssic"}, \code{"aicc"}, \code{"hybrid"}.
#' @param gamma exponent for the sSIC penalty.
#' @return A list with `s` (criterion value), `llh`, `p` (penalty), `pv`
#'   (speed-penalty component), `logical`.
#' @export
CS <- function(t, x, y, r, p1, s_cap, with_speed, sd = NA, pen = "ssic", gamma) {
  n <- length(t)
  cps <- which(r == 1) + 1
  pla <- piecewise_linear_con(t, x, y, cps)

  if (!isTRUE(pla$logical)) {
    return(list(logical = FALSE))
  }

  llh <- loglikelihood(n, pla$path_RSS, sd = sd)
  seg_vel <- pla$path_segvel
  pv <- sum(pmax(0, seg_vel - s_cap))
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
