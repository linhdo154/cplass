#' @keywords internal
#' @noRd
d_fn <- function(k, p) {
  as.integer(k > p)
}

# ---- Type 1: independent changepoint vector (q_new) ------------------------

#' @rdname proposals
#' @export
log_q_new_pmf <- function(lambda_r, dt, N, r) {
  K_r <- sum(r)
  K_r * log(1 - exp(-lambda_r * dt)) - (N - 2 - K_r) * (lambda_r * dt)
}

#' @rdname proposals
#' @export
q_new <- function(lambda_r, dt, N) {
  p <- 1 - exp(-lambda_r * dt)
  stats::rbinom(N - 2, size = 1, prob = p)
}

# ---- Type 2: birth/death of a single changepoint (q_bd) --------------------

#' @rdname proposals
#' @export
log_q_bd_pmf <- function(r, status) {
  N_r <- length(r)
  K_r <- sum(r)
  if (status == 0) {
    if (!any(r == 1)) {
      return(log(0))
    }
    return(log(1 / (2 * K_r)))
  }
  if (status == 1) {
    if (!any(r == 0)) {
      return(log(0))
    }
    return(log(1 / (2 * (N_r - K_r))))
  }
  stop("`status` must be 0 or 1", call. = FALSE)
}

#' @rdname proposals
#' @export
q_bd <- function(r_cur) {
  r_prop <- r_cur
  N_r <- length(r_cur)

  if (!any(r_cur == 0) || !any(r_cur == 1)) {
    s <- sample.int(N_r, size = 1)
    r_prop[s] <- 1 - r_cur[s]
    status <- r_prop[s]
  } else {
    status <- sample(c(0, 1), size = 1)
    if (status == 1) {
      s <- sample(which(r_cur == 0), size = 1)
      r_prop[s] <- 1
    } else {
      s <- sample(which(r_cur == 1), size = 1)
      r_prop[s] <- 0
    }
  }

  list(r_prop = r_prop, status = status)
}

# ---- Type 4: single changepoint location shift (q_shift) -------------------

#' @rdname proposals
#' @export
q_shift <- function(r_cur) {
  ones_idx <- which(r_cur == 1)
  zeros_idx <- which(r_cur == 0)

  if (length(ones_idx) == 0 || length(zeros_idx) == 0) {
    stop("Need at least one 1 and one 0 in the changepoint vector.", call. = FALSE)
  }

  s <- if (length(ones_idx) == 1) ones_idx else sample(ones_idx, size = 1)
  ss <- if (length(zeros_idx) == 1) zeros_idx else sample(zeros_idx, size = 1)

  r_prop <- r_cur
  r_prop[c(s, ss)] <- 1 - r_cur[c(s, ss)]
  r_prop
}

# ---- Type 3: birth/death of a segment (paired changepoints) ---------------

#' @rdname proposals
#' @export
q_as_pmf <- function(r) {
  N_r <- length(r)
  K_r <- sum(r)
  cp <- sort(c(1, which(r == 1) + 1, N_r + 1))
  d <- diff(cp)
  com <- ifelse(d < 2, 0, (d - 1) * (d - 2))
  (1 / (2 * ((N_r - K_r) * (N_r - K_r - 1)))) * sum(com)
}

#' @rdname proposals
#' @export
q_as <- function(r_cur) {
  r_prop <- r_cur
  N_r <- length(r_cur)
  K_r <- sum(r_cur)
  idx_add <- sample.int(K_r + 1, size = 1)
  cp <- sort(c(1, which(r_cur == 1) + 1, N_r + 2))

  if (cp[idx_add + 1] - cp[idx_add] > 2) {
    new <- sample(cp[idx_add]:cp[idx_add + 1], 2)
    cp_prop <- unique(sort(c(which(r_cur == 1) + 1, new)))
    r_prop[cp_prop - 1] <- 1
  }
  r_prop
}

#' @rdname proposals
#' @export
q_ds_pmf <- function(r) {
  K_r <- sum(r)
  if (K_r >= 2) 1 / K_r else 0
}

#' @rdname proposals
#' @export
q_ds <- function(r_cur) {
  r_prop <- r_cur
  ones <- which(r_cur == 1)

  if (r_cur[length(r_cur)] == 1 && length(ones) <= 2) {
    s <- ones[-length(ones)]
  } else {
    s <- sample(ones[-length(ones)], size = 1)
  }

  r_prop[s] <- 0
  r_prop[ones[which(ones == s) + 1]] <- 0
  r_prop
}

#' @rdname proposals
#' @export
q_bd2 <- function(r_cur) {
  if (!any(r_cur == 1)) {
    s <- sample.int(length(r_cur), 2)
    r_prop <- r_cur
    r_prop[s] <- 1
    status <- 1
  } else if (!any(r_cur == 0)) {
    s <- sample.int(length(r_cur), 1)
    r_prop <- r_cur
    r_prop[s] <- 0
    if (s == length(r_cur)) r_prop[s - 1] <- 0 else r_prop[s + 1] <- 0
    status <- 0
  } else if (length(which(r_cur == 1)) == 1) {
    r_prop <- q_as(r_cur)
    status <- 1
  } else {
    status <- sample(c(0, 1), 1)
    r_prop <- if (status == 1) q_as(r_cur) else q_ds(r_cur)
  }

  list(r_prop = r_prop, status = status)
}

#' @rdname proposals
#' @export
log_q_bd2_pmf <- function(r, status) {
  p <- if (status == 0) q_ds_pmf(r) else q_as_pmf(r)
  log(p)
}

# ---- Combined proposal ------------------------------------------------------

#' Proposal functions for the CPLASS Metropolis-Hastings sampler
#'
#' These implement the four proposal types described in Section 2.4 of the
#' companion paper: (1) an independent changepoint draw (`q_new`), (2)
#' birth/death of a single changepoint (`q_bd`), (3) birth/death of a
#' segment / paired changepoints (`q_bd2`, via `q_as`/`q_ds`), and (4) a
#' location shift of a single changepoint (`q_shift`). `proposal_function()`
#' combines them with mixture weights `u1, u2, u3`; `pproposal()` evaluates
#' the corresponding proposal density for the Metropolis-Hastings ratio.
#'
#' @param u draw from Uniform(0,1) selecting which proposal type to use.
#' @param r_cur,r current changepoint indicator vector.
#' @param N number of observations.
#' @param lambda_r rate parameter for the Type-1 (independent) proposal.
#' @param dt observation time step.
#' @name proposals
#' @export
proposal_function <- function(u, r_cur, N, lambda_r, dt) {
  saturated <- !any(r_cur == 0) || !any(r_cur == 1)

  if (saturated) {
    if (u <= 1 / 4) {
      return(list(r_prop = q_new(lambda_r, dt, N), status = 2))
    } else if (u <= 1 / 2) {
      pp <- q_bd(r_cur)
      return(list(r_prop = pp$r_prop, status = pp$status))
    } else {
      pp <- q_bd2(r_cur)
      return(list(r_prop = pp$r_prop, status = pp$status))
    }
  }

  if (u <= 1 / 4) {
    list(r_prop = q_new(lambda_r, dt, N), status = 2)
  } else if (u <= 3 / 8) {
    pp <- q_bd(r_cur)
    list(r_prop = pp$r_prop, status = pp$status)
  } else if (u <= 1 / 2) {
    pp <- q_bd2(r_cur)
    list(r_prop = pp$r_prop, status = pp$status)
  } else {
    list(r_prop = q_shift(r_cur), status = 3)
  }
}

#' @rdname proposals
#' @param r_given the changepoint vector whose proposal density is being
#'   evaluated (either the proposed or current vector, depending on
#'   direction of the MH ratio).
#' @param status proposal-type-specific status flag returned alongside
#'   `r_prop` by `q_bd()`/`q_bd2()`.
#' @export
pproposal <- function(u, r, N, lambda_r, r_given, status, dt) {
  saturated <- !any(r_given == 0) || !any(r_given == 1)

  if (saturated) {
    if (u <= 1 / 4) {
      log_q_new_pmf(lambda_r, dt, N, r)
    } else if (u <= 1 / 2) {
      log_q_bd_pmf(r_given, status)
    } else {
      log_q_bd2_pmf(r_given, status)
    }
  } else {
    if (u <= 1 / 4) {
      log_q_new_pmf(lambda_r, dt, N, r)
    } else if (u <= 3 / 8) {
      log_q_bd_pmf(r_given, status)
    } else if (u <= 1 / 2) {
      log_q_bd2_pmf(r_given, status)
    } else {
      log(1 / (N - 2))
    }
  }
}
