test_that("CPLASS() runs end-to-end on a synthetic 2-changepoint trajectory", {
  set.seed(42)
  n <- 150
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  seg <- findInterval(seq_along(t), c(0, 50, 100, n + 1))
  v <- list(c(0.05, -0.02), c(0, 0.25), c(-0.05, 0.02))
  x <- y <- numeric(n)
  for (i in 2:n) {
    vv <- v[[seg[i]]]
    x[i] <- x[i - 1] + vv[1] * 0.05
    y[i] <- y[i - 1] + vv[2] * 0.05
  }
  x <- x + rnorm(n, 0, 0.01)
  y <- y + rnorm(n, 0, 0.01)

  res <- CPLASS(t, x, y, iter_max = 300, burn_in = 50, show_progress = FALSE)

  expect_type(res, "list")
  expect_true(all(c("segments_inferred", "path_inferred") %in% names(res)))
  expect_true(nrow(res$segments_inferred) >= 1)
  expect_equal(nrow(res$path_inferred), n)
  expect_true(all(res$segments_inferred$speeds >= 0))
})

test_that("CPLASS() errors informatively when burn_in >= iter_max", {
  set.seed(1)
  n <- 60
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  x <- cumsum(rnorm(n, 0, 0.02))
  y <- cumsum(rnorm(n, 0, 0.02))

  expect_message(
    CPLASS(t, x, y, iter_max = 50, burn_in = 50, show_progress = FALSE),
    "failed too many times|failed, retrying"
  )
})

test_that("MHsearch acceptance rate is between 0 and 1 and info_table has the right shape", {
  set.seed(5)
  n <- 100
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  x <- cumsum(rnorm(n, 0, 0.02))
  y <- cumsum(rnorm(n, 0, 0.02))
  dt <- min(diff(t))

  out <- MHsearch(t, x, y, dt, iter_max = 200, burn_in = 50, show_progress = FALSE)
  expect_equal(nrow(out$info_table), 150)
  expect_length(out$cps_list, 150)
  expect_true(out$update_info$accept_rate >= 0 && out$update_info$accept_rate <= 1)
})
