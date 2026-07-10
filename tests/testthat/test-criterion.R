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
