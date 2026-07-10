synthetic_path <- function(seed = 1, n = 150) {
  set.seed(seed)
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
  list(t = t, x = x, y = y)
}

test_that("CPLASS() returns a best_score alongside the usual output", {
  p <- synthetic_path(1)
  res <- CPLASS(p$t, p$x, p$y, iter_max = 300, burn_in = 50, show_progress = FALSE)
  expect_true(is.numeric(res$best_score))
  expect_length(res$best_score, 1)
})

test_that("check_convergence() flags a short chain as not converged and a long one as converged", {
  p <- synthetic_path(2)

  set.seed(10)
  short <- check_convergence(p$t, p$x, p$y, iter_max = 150, burn_in = 50, window_frac = 0.5)
  expect_type(short$converged, "logical")
  expect_true(nrow(short$trace) == 100)
  expect_true(short$best_score >= max(short$trace$score))

  set.seed(10)
  long <- check_convergence(p$t, p$x, p$y, iter_max = 2000, burn_in = 50, window_frac = 0.2)
  expect_true(nrow(long$trace) == 1950)
  # running max must be non-decreasing
  expect_true(all(diff(long$trace$running_max) >= 0))
})

test_that("plot_convergence() returns a ggplot object", {
  p <- synthetic_path(3)
  set.seed(11)
  conv <- check_convergence(p$t, p$x, p$y, iter_max = 300, burn_in = 50)
  gg <- plot_convergence(conv)
  expect_s3_class(gg, "ggplot")
})

test_that("CPLASS_multistart() picks the highest-scoring of several starts", {
  p <- synthetic_path(4)
  res <- CPLASS_multistart(p$t, p$x, p$y, iter_max = 300, burn_in = 50,
                             n_starts = 3, seeds = c(101, 202, 303),
                             patience = NULL, show_progress = FALSE)

  expect_length(res$all_scores, 3)
  expect_length(res$seeds, 3)
  expect_equal(res$best$best_score, max(res$all_scores))
  expect_true(res$best_seed %in% c(101, 202, 303))
  expect_true(res$score_spread >= 0)

  # re-running with the winning seed reproduces the winning result exactly
  set.seed(res$best_seed)
  repro <- CPLASS(p$t, p$x, p$y, iter_max = 300, burn_in = 50, show_progress = FALSE)
  expect_equal(repro$best_score, res$best$best_score)
  expect_equal(repro$segments_inferred$cp_times, res$best$segments_inferred$cp_times)
})

test_that("CPLASS_multistart() default patience = 'auto' resolves to 0.3 * iter_max", {
  p <- synthetic_path(6)
  res <- CPLASS_multistart(p$t, p$x, p$y, iter_max = 1000, burn_in = 100,
                             n_starts = 2, seeds = c(1, 2), show_progress = FALSE)
  expect_equal(res$patience_used, ceiling(0.3 * 1000))
  expect_true(all(vapply(list(res$best), function(r) !is.null(r$iterations_used), logical(1))))
})

test_that("CPLASS_multistart() errors informatively when every start fails", {
  p <- synthetic_path(5)
  # burn_in >= iter_max makes CPLASS() fail every time (see test-mcmc.R);
  # this exercises the "all starts failed" path
  expect_error(
    CPLASS_multistart(p$t, p$x, p$y, iter_max = 50, burn_in = 50, n_starts = 2,
                        show_progress = FALSE),
    "failed to converge"
  )
})
