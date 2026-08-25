test_that("case 1: c(r) = n-2 (saturated) is rejected", {
  n <- 6
  r <- rep(1L, n - 2) # c(r) = 4 = n - 2
  expect_false(.r_is_admissible(r, n))
})

test_that("case 2: c(r) = n-3 is structurally allowed (all else satisfied)", {
  n <- 6
  r <- c(1L, 1L, 1L, 0L) # c(r) = 3 = n - 3, K(r) = 4 = n - 2 (at the ceiling)
  expect_true(.r_is_admissible(r, n))
})

test_that("case 3: exceeding a supplied K_max is rejected", {
  n <- 30
  r <- integer(n - 2)
  r[c(5, 15)] <- 1L # c(r) = 2, K(r) = 3
  expect_true(.r_is_admissible(r, n)) # allowed with no ceiling supplied
  expect_false(.r_is_admissible(r, n, K_max = 2)) # K(r) = 3 > K_max = 2
  expect_true(.r_is_admissible(r, n, K_max = 3)) # K(r) = 3 <= K_max = 3
})

test_that("case 4: an ordinary well-separated configuration is accepted", {
  n <- 100
  r <- integer(n - 2)
  r[c(10, 40, 70)] <- 1L
  expect_true(.r_is_admissible(r, n))
  expect_true(.r_is_admissible(r, n, K_max = 10, min_gap = 1))
})

test_that("case 5: rank-deficient design is still rejected by the existing QR check", {
  # Regression check: confirms the pre-existing full-rank guard in
  # .fit_pla() is unchanged by this patch. Duplicate/adjacent "changepoint"
  # indices make the hinge design rank-deficient.
  t <- seq(0, 1, length.out = 10)
  set.seed(1)
  x <- rnorm(10)
  y <- rnorm(10)
  fit <- piecewise_linear_con(t, x, y, c(5L, 5L))
  expect_false(isTRUE(fit$logical))
})

test_that("case 6: nonpositive residual degrees of freedom is rejected safely", {
  # n = 10, cp = 2:9 (8 changepoints) makes the hinge design SQUARE:
  # ncol(A) = 2 + length(cp) = 10 = n. This design is full column rank (so
  # the pre-existing QR check alone would NOT reject it) but interpolates
  # the data exactly, giving RSS = 0 and resid_df = 2*10 - 2*9 - 2 = 0.
  # This is exactly the unbounded-criterion failure mode Comment 2
  # describes, and is the case the new check exists to catch.
  n <- 10
  t <- seq(0, (n - 1) * 0.05, length.out = n)
  set.seed(2)
  x <- cumsum(rnorm(n, 0, 0.02))
  y <- cumsum(rnorm(n, 0, 0.02))
  cp <- 2:9
  expect_equal(length(cp), n - 2) # c(r) = n - 2, saturated

  # Confirm the design really is full rank, i.e. that the pre-existing QR
  # check alone would have let this through.
  hinge <- sapply(cp, function(c) pmax(t - t[c], 0))
  A <- cbind(1, t, hinge)
  expect_equal(qr(A)$rank, ncol(A))

  fit <- piecewise_linear_con(t, x, y, cp)
  expect_false(isTRUE(fit$logical))

  # The upstream admissibility guard should also reject this r (defense in
  # depth: this configuration should never reach .fit_pla via CS() at all).
  r <- cp_to_r(cp, n)
  expect_false(.r_is_admissible(r, n))
})

test_that("case 7: (numerically) zero RSS is rejected safely", {
  # Same exact-interpolation construction as case 6; RSS = 0 and
  # resid_df = 0 co-occur structurally at this boundary (n = ncol(A)
  # exactly always gives an exact interpolating fit), so this is verified
  # via the same scenario rather than an artificially separate one.
  n <- 8
  t <- seq(0, (n - 1) * 0.05, length.out = n)
  set.seed(3)
  x <- cumsum(rnorm(n, 0, 0.02))
  y <- cumsum(rnorm(n, 0, 0.02))
  cp <- 2:7 # length = n - 2 = 6, ncol(A) = 8 = n
  fit <- piecewise_linear_con(t, x, y, cp)
  expect_false(isTRUE(fit$logical))
})

test_that("case 8: no NaN/Inf sigma_hat or criterion leaks from an invalid configuration", {
  n <- 10
  t <- seq(0, (n - 1) * 0.05, length.out = n)
  set.seed(4)
  x <- cumsum(rnorm(n, 0, 0.02))
  y <- cumsum(rnorm(n, 0, 0.02))
  cp <- 2:9

  fit <- piecewise_linear_con(t, x, y, cp)
  expect_false(isTRUE(fit$logical))
  # An invalid fit returns only `logical = FALSE`; no sigma_hat/criterion
  # value is present to leak a NaN/Inf through.
  expect_null(fit$sigma_hat)
  expect_null(fit$path_segeta)

  # Same check one level up, through CS()/the admissibility guard.
  r <- cp_to_r(cp, n)
  res <- CS(t, x, y, r, p1 = 1, s_cap = 1, with_speed = 0, gamma = 1.01)
  expect_false(isTRUE(res$logical))
  expect_null(res$s)
})

test_that("malformed r is rejected without erroring", {
  n <- 20
  expect_false(.r_is_admissible(rep(1, n - 3), n)) # wrong length
  expect_false(.r_is_admissible(c(0, 1, 2, rep(0, n - 5)), n)) # non-binary entry
  expect_false(.r_is_admissible(c(NA, rep(0, n - 3)), n)) # NA entry
})

test_that("min_gap = 1 (default) is non-binding: adjacent changepoints allowed", {
  n <- 20
  r <- integer(n - 2)
  r[c(5, 6)] <- 1L # adjacent positions -> adjacent changepoints
  expect_true(.r_is_admissible(r, n)) # min_gap = 1 default: not rejected
  expect_true(.r_is_admissible(r, n, min_gap = 1))
})

test_that("min_gap = 2 rejects an adjacent-changepoint configuration", {
  n <- 20
  r <- integer(n - 2)
  r[c(5, 6)] <- 1L # adjacent positions -> gap of 1 between them
  expect_true(.r_is_admissible(r, n, min_gap = 1)) # not rejected at default
  expect_false(.r_is_admissible(r, n, min_gap = 2)) # gap of 1 < min_gap = 2

  # a configuration with all gaps >= 2 should still be accepted under min_gap = 2
  r2 <- integer(n - 2)
  r2[c(5, 8, 12)] <- 1L
  expect_true(.r_is_admissible(r2, n, min_gap = 2))
})

test_that("K_max effective ceiling is min(K_max, n-2), not K_max alone", {
  n <- 8 # n - 2 = 6, so K(r) can be at most 6 structurally (c(r) <= n-3 = 5)
  r <- integer(n - 2)
  r[1:5] <- 1L # c(r) = 5 = n - 3, K(r) = 6 = n - 2 (at the structural ceiling)
  expect_true(.r_is_admissible(r, n)) # no K_max supplied
  # A K_max looser than the structural ceiling changes nothing
  expect_true(.r_is_admissible(r, n, K_max = 100))
  # A K_max tighter than the structural ceiling does bind
  expect_false(.r_is_admissible(r, n, K_max = 5))
})

test_that("K_max/min_gap parameter validation rejects invalid values with a clear error", {
  N <- 20
  t <- seq(0, (N - 1) * 0.05, by = 0.05)
  x <- cumsum(rnorm(N, 0, 0.02))
  y <- cumsum(rnorm(N, 0, 0.02))

  expect_error(MHsearch(t, x, y, dt = 0.05, iter_max = 10, burn_in = 2, K_max = 0), "K_max")
  expect_error(MHsearch(t, x, y, dt = 0.05, iter_max = 10, burn_in = 2, K_max = -1), "K_max")
  expect_error(MHsearch(t, x, y, dt = 0.05, iter_max = 10, burn_in = 2, K_max = c(4, 5)), "K_max")
  expect_error(MHsearch(t, x, y, dt = 0.05, iter_max = 10, burn_in = 2, K_max = NA), "K_max")
  expect_error(MHsearch(t, x, y, dt = 0.05, iter_max = 10, burn_in = 2, min_gap = 0), "min_gap")
  expect_error(MHsearch(t, x, y, dt = 0.05, iter_max = 10, burn_in = 2, min_gap = -2), "min_gap")
  expect_error(MHsearch(t, x, y, dt = 0.05, iter_max = 10, burn_in = 2, min_gap = NA), "min_gap")
  expect_error(MHsearch(t, x, y, dt = 0.05, iter_max = 10, burn_in = 2, min_gap = c(1, 2)), "min_gap")

  # valid values should NOT error at the validation step
  expect_true(.validate_K_max_min_gap(NULL, 1L))
  expect_true(.validate_K_max_min_gap(5, 2))
})
