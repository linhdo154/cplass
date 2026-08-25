test_that("sSIC and AICc are monotonic in the number of changepoints", {
  n <- 500
  ncp_vals <- 0:10
  ssic_vals <- sapply(ncp_vals, sSIC, n = n, gamma = 1.01)
  aicc_vals <- sapply(ncp_vals, AICc, n = n)
  expect_true(all(diff(ssic_vals) > 0))
  expect_true(all(diff(aicc_vals) > 0))
})

test_that("CS() returns a lower (or equal) criterion for a subset of changepoints, up to penalty", {
  set.seed(3)
  n <- 100
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  x <- cumsum(rnorm(n, 0, 0.02))
  y <- cumsum(rnorm(n, 0, 0.02))

  r0 <- integer(n - 2)
  r1 <- r0
  r1[40] <- 1

  s0 <- CS(t, x, y, r0, p1 = 1, s_cap = 1, with_speed = 0, gamma = 1.01)
  s1 <- CS(t, x, y, r1, p1 = 1, s_cap = 1, with_speed = 0, gamma = 1.01)

  expect_true(s0$logical)
  expect_true(s1$logical)
  # adding a changepoint should always improve (or match) the raw log-likelihood
  expect_true(s1$llh >= s0$llh - 1e-8)
})

test_that("cp_to_r and r_to_cp round-trip", {
  n <- 50
  cp <- c(5L, 12L, 30L)
  r <- cp_to_r(cp, n)
  expect_equal(sort(r_to_cp(r)), sort(cp))
})

test_that("relative-excess speed penalty is invariant to spatial units", {
  set.seed(41)
  n <- 81L
  t <- seq(0, 4, length.out = n)
  x <- ifelse(t <= 2, 0, 4 * (t - 2)) + rnorm(n, sd = 0.001)
  y <- rnorm(n, sd = 0.001)
  r <- cp_to_r(41L, n)

  um <- CS(t, x, y, r, p1 = 1, s_cap = 2, with_speed = 1,
           eta = 1, gamma = 1.01)
  nm <- CS(t, 1000 * x, 1000 * y, r, p1 = 1, s_cap = 2000,
           with_speed = 1, eta = 1, gamma = 1.01)

  expect_true(um$logical)
  expect_true(nm$logical)
  expect_gt(um$pv, 0)
  expect_equal(nm$pv, um$pv, tolerance = 1e-8)
})

test_that("eta controls an extreme short high-speed segment", {
  set.seed(42)
  n <- 101L
  t <- seq(0, 5, length.out = n)
  anchor <- 20 * pmax(0, pmin(t - 2.45, 0.10))
  x <- anchor + rnorm(n, sd = 0.005)
  y <- rnorm(n, sd = 0.005)
  r <- cp_to_r(c(50L, 52L), n)

  no_speed <- CS(t, x, y, r, p1 = 1, s_cap = 2, with_speed = 1,
                 eta = 0, gamma = 1.01)
  regularized <- CS(t, x, y, r, p1 = 1, s_cap = 2, with_speed = 1,
                    eta = 5, gamma = 1.01)

  expect_true(no_speed$logical)
  expect_true(regularized$logical)
  expect_gt(regularized$pv, 0)
  expect_lt(regularized$s, no_speed$s)
  expect_equal(no_speed$s - regularized$s, regularized$pv,
               tolerance = 1e-8)
})

test_that("speed-penalty inputs are validated", {
  n <- 20L
  t <- seq_len(n)
  r <- integer(n - 2L)

  expect_error(
    CS(t, t, t, r, p1 = 1, s_cap = 0, with_speed = 1,
       eta = 1, gamma = 1.01),
    "positive"
  )
  expect_error(
    CS(t, t, t, r, p1 = 1, s_cap = 1, with_speed = 1,
       eta = -1, gamma = 1.01),
    "non-negative"
  )
})
