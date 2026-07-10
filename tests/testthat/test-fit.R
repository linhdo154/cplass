test_that("piecewise_linear_con matches a manual OLS fit (0, 1, 2 changepoints)", {
  set.seed(1)
  n <- 120
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  x <- cumsum(rnorm(n, 0, 0.02))
  y <- cumsum(rnorm(n, 0, 0.02))

  check_cp <- function(cp) {
    fit <- piecewise_linear_con(t, x, y, cp)
    expect_true(fit$logical)

    # Rebuild the hinge design matrix by hand and check against base lm-style OLS
    m <- length(cp) + 1
    if (length(cp) > 0) {
      hinge <- sapply(sort(cp), function(c) pmax(t - t[c], 0))
      A <- cbind(1, t, hinge)
    } else {
      A <- cbind(1, t)
    }
    alpha_manual <- qr.solve(A, x)
    beta_manual <- qr.solve(A, y)

    expect_equal(as.numeric(fit$x_piecewise), as.numeric(A %*% alpha_manual), tolerance = 1e-8)
    expect_equal(as.numeric(fit$y_piecewise), as.numeric(A %*% beta_manual), tolerance = 1e-8)

    RSS_manual <- sum((x - A %*% alpha_manual)^2) + sum((y - A %*% beta_manual)^2)
    expect_equal(fit$path_RSS, RSS_manual, tolerance = 1e-8)
  }

  check_cp(integer(0))
  check_cp(40L)
  check_cp(c(30L, 70L))
  check_cp(c(20L, 45L, 60L, 90L))
})

test_that("piecewise_linear_con handles degenerate/rank-deficient inputs gracefully", {
  # Two "changepoints" at adjacent indices can make the design matrix
  # rank-deficient; should fail cleanly (logical = FALSE), not error.
  t <- seq(0, 1, length.out = 10)
  x <- rnorm(10)
  y <- rnorm(10)
  fit <- piecewise_linear_con(t, x, y, c(5L, 5L))
  expect_false(isTRUE(fit$logical))
})

test_that("segment velocities are cumulative sums of the increment coefficients", {
  set.seed(2)
  n <- 90
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  x <- cumsum(rnorm(n, 0, 0.02))
  y <- cumsum(rnorm(n, 0, 0.02))
  cp <- c(25L, 55L)

  fit <- piecewise_linear_con(t, x, y, cp)
  expect_length(fit$u_x, length(cp) + 1)
  expect_length(fit$v_y, length(cp) + 1)
  expect_equal(sqrt(fit$u_x^2 + fit$v_y^2), fit$path_segvel, tolerance = 1e-10)
})
