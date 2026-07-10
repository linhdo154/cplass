rbern <- function(n, prob) rbinom(n, size = 1, prob = prob)
orig_env <- new.env()
local({
  real_library <- base::library
  library <- function(pkg, ...) {
    pkg_name <- as.character(substitute(pkg))
    if (requireNamespace(pkg_name, quietly = TRUE)) real_library(pkg_name, character.only = TRUE)
  }
  assign("library", library, envir = orig_env)
  assign("rbern", rbern, envir = orig_env)
  sys.source("/home/claude/CPLASS/code/PENALTY.R", envir = orig_env)
  sys.source("/home/claude/CPLASS/code/CPLASS.R", envir = orig_env)
}, envir = orig_env)
orig_pla <- get("piecewise_linear_con", envir = orig_env)

suppressMessages(library(cplass))

set.seed(1)
for (n in c(100, 500, 2000)) {
  t <- seq(0, (n-1)*0.05, by=0.05)
  x <- cumsum(rnorm(n, 0, 0.05))
  y <- cumsum(rnorm(n, 0, 0.05))
  cps <- sort(sample(3:(n-3), min(8, floor(n/20))))

  reps <- 200
  t_old <- system.time(for (k in 1:reps) orig_pla(t, x, y, cps))["elapsed"]
  t_new <- system.time(for (k in 1:reps) cplass::piecewise_linear_con(t, x, y, cps))["elapsed"]

  cat(sprintf("n=%5d | old: %7.4fs (%d calls) | new: %7.4fs | speedup: %.1fx\n",
              n, t_old, reps, t_new, t_old/t_new))
}
