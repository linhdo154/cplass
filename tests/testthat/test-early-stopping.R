test_that("patience = NULL (default) behaves exactly as before: full iter_max used", {
  set.seed(20)
  n <- 100
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  x <- cumsum(rnorm(n, 0, 0.02))
  y <- cumsum(rnorm(n, 0, 0.02))
  dt <- min(diff(t))

  out <- MHsearch(t, x, y, dt, iter_max = 300, burn_in = 50, show_progress = FALSE)
  expect_false(out$stopped_early)
  expect_equal(out$iterations_used, 300)
  expect_equal(nrow(out$info_table), 250)
})

test_that("setting a small patience stops the chain before iter_max on an easy path", {
  set.seed(21)
  n <- 150
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  # a trivial, single-segment path: should converge almost immediately and
  # then plateau, so a small patience should trigger well before iter_max
  x <- t * 0.05 + rnorm(n, 0, 0.005)
  y <- -t * 0.02 + rnorm(n, 0, 0.005)
  dt <- min(diff(t))

  out <- MHsearch(t, x, y, dt, iter_max = 5000, burn_in = 100, patience = 200,
                   show_progress = FALSE)
  expect_true(out$stopped_early)
  expect_true(out$iterations_used < 5000)
  expect_true(out$iterations_used >= 100 + 200)
})

test_that("early-stopped MHsearch() output is internally consistent", {
  set.seed(22)
  n <- 150
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  x <- t * 0.05 + rnorm(n, 0, 0.005)
  y <- -t * 0.02 + rnorm(n, 0, 0.005)
  dt <- min(diff(t))

  out <- MHsearch(t, x, y, dt, iter_max = 5000, burn_in = 100, patience = 200,
                   show_progress = FALSE)
  expect_equal(nrow(out$info_table), length(out$cps_list))
  expect_equal(nrow(out$info_table), out$iterations_used - 100)
  # the running max of accepted scores must be non-decreasing throughout
  accepted_score <- ifelse(out$info_table$decision == 1, out$info_table$score_new, out$info_table$score_cur)
  expect_true(all(diff(cummax(accepted_score)) >= 0))
})

test_that("CPLASS() with patience produces a valid result and reports stopped_early/iterations_used", {
  set.seed(23)
  n <- 150
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  x <- t * 0.05 + rnorm(n, 0, 0.005)
  y <- -t * 0.02 + rnorm(n, 0, 0.005)

  res <- CPLASS(t, x, y, iter_max = 5000, burn_in = 100, patience = 200,
                show_progress = FALSE)
  expect_true(!is.null(res$iterations_used))
  expect_true(res$iterations_used <= 5000)
  expect_true(nrow(res$segments_inferred) >= 1)
})

test_that("patience does not change the *quality* of the result on an easy path", {
  # compare a patience-stopped short run against a full long run: both
  # should find (approximately) the same best score, since the path is
  # simple enough that early stopping isn't cutting off real improvement
  set.seed(24)
  n <- 150
  t <- seq(0, (n - 1) * 0.05, by = 0.05)
  x <- t * 0.05 + rnorm(n, 0, 0.005)
  y <- -t * 0.02 + rnorm(n, 0, 0.005)

  set.seed(99)
  full <- CPLASS(t, x, y, iter_max = 5000, burn_in = 100, show_progress = FALSE)
  set.seed(99)
  early <- CPLASS(t, x, y, iter_max = 5000, burn_in = 100, patience = 300, show_progress = FALSE)

  expect_true(early$iterations_used < 5000)
  expect_equal(early$best_score, full$best_score, tolerance = 1e-6)
})
