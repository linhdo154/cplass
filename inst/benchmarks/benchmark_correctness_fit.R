rbern <- function(n, prob) rbinom(n, size = 1, prob = prob)
orig_env <- new.env()
local({
  real_library <- base::library
  library <- function(pkg, ...) {
    pkg_name <- as.character(substitute(pkg))
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      real_library(pkg_name, character.only = TRUE)
    }
  }
  assign("library", library, envir = orig_env)
  assign("rbern", rbern, envir = orig_env)
  sys.source("/home/claude/CPLASS/code/PENALTY.R", envir = orig_env)
  sys.source("/home/claude/CPLASS/code/CPLASS.R", envir = orig_env)
}, envir = orig_env)

orig_pla <- get("piecewise_linear_con", envir = orig_env)

suppressMessages(library(cplass))

set.seed(42)
n <- 300
t <- seq(0, (n - 1) * 0.05, by = 0.05)
v <- list(c(0.1, -0.05), c(0, 0.2), c(-0.1, 0.05))
cp_true <- c(100, 200)
seg <- findInterval(seq_along(t), c(0, cp_true, n + 1))
x <- y <- numeric(n)
for (i in 2:n) {
  vv <- v[[seg[i]]]
  x[i] <- x[i - 1] + vv[1] * 0.05
  y[i] <- y[i - 1] + vv[2] * 0.05
}
x <- x + rnorm(n, 0, 0.02)
y <- y + rnorm(n, 0, 0.02)

cps_to_test <- list(integer(0), 50L, c(100L, 200L), c(30L, 90L, 180L, 250L))

max_abs_diff <- function(a, b) max(abs(a - b))

all_ok <- TRUE
for (cps in cps_to_test) {
  o <- orig_pla(t, x, y, cps)
  n_ <- cplass::piecewise_linear_con(t, x, y, cps)

  diffs <- c(
    RSS         = max_abs_diff(o$path_RSS, n_$path_RSS),
    segvel      = max_abs_diff(o$path_segvel, n_$path_segvel),
    segtime     = max_abs_diff(o$path_segtime, n_$path_segtime),
    segtheta    = max_abs_diff(o$path_segtheta, n_$path_segtheta),
    segeta      = max_abs_diff(o$path_segeta, n_$path_segeta),
    u_x         = max_abs_diff(o$u_x, n_$u_x),
    v_y         = max_abs_diff(o$v_y, n_$v_y),
    x_piecewise = max_abs_diff(o$x_piecewise, n_$x_piecewise),
    y_piecewise = max_abs_diff(o$y_piecewise, n_$y_piecewise)
  )
  ok <- all(diffs < 1e-8)
  all_ok <- all_ok && ok
  cat("cps =", paste(cps, collapse = ","), " -> max abs diffs:\n")
  print(diffs)
  cat("MATCH:", ok, "\n\n")
}

cat("=== ALL MATCH:", all_ok, "===\n")
