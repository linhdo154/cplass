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
orig_CPLASS <- get("CPLASS", envir = orig_env)

suppressMessages(library(cplass))
suppressMessages(library(readr))
suppressMessages(library(dplyr))

data <- suppressMessages(read_csv("/home/claude/CPLASS/data/Real_21_Periphery.csv", show_col_types = FALSE))
path1 <- data %>% filter(index_path == 1)
cat("=== REAL DATA: path1, n obs =", nrow(path1), "===\n\n")

for (iter_max in c(1000, 3000)) {
  set.seed(123)
  tm_old <- system.time(res_old <- orig_CPLASS(path1$t, path1$x, path1$y, iter_max = iter_max, burn_in = 100))["elapsed"]
  cat(sprintf("%-8s | iter_max=%5d | time=%7.3fs | n_segments=%d\n", "OLD", iter_max, tm_old, nrow(res_old$segments_inferred)))

  set.seed(123)
  tm_new <- system.time(res_new <- cplass::CPLASS(path1$t, path1$x, path1$y, iter_max = iter_max, burn_in = 100, show_progress = FALSE))["elapsed"]
  cat(sprintf("%-8s | iter_max=%5d | time=%7.3fs | n_segments=%d\n", "NEW", iter_max, tm_new, nrow(res_new$segments_inferred)))
  cat(sprintf("  -> speedup: %.1fx\n\n", tm_old / tm_new))
}

# --- synthetic data, larger n, matching paper's simulation style ---
cat("=== SYNTHETIC DATA (n=600, 2 true changepoints, 20Hz) ===\n\n")
set.seed(7)
n <- 600
t <- seq(0, (n-1)*0.05, by = 0.05)
seg <- findInterval(seq_along(t), c(0, 200, 400, n+1))
v <- list(c(0.05,-0.02), c(0,0.25), c(-0.05,0.02))
x <- y <- numeric(n)
for (i in 2:n) { vv <- v[[seg[i]]]; x[i] <- x[i-1]+vv[1]*0.05; y[i] <- y[i-1]+vv[2]*0.05 }
x <- x + rnorm(n, 0, 0.02); y <- y + rnorm(n, 0, 0.02)

for (iter_max in c(1000, 3000)) {
  set.seed(99)
  tm_old <- system.time(res_old <- orig_CPLASS(t, x, y, iter_max = iter_max, burn_in = 100))["elapsed"]
  cat(sprintf("%-8s | iter_max=%5d | time=%7.3fs | n_segments=%d\n", "OLD", iter_max, tm_old, nrow(res_old$segments_inferred)))

  set.seed(99)
  tm_new <- system.time(res_new <- cplass::CPLASS(t, x, y, iter_max = iter_max, burn_in = 100, show_progress = FALSE))["elapsed"]
  cat(sprintf("%-8s | iter_max=%5d | time=%7.3fs | n_segments=%d\n", "NEW", iter_max, tm_new, nrow(res_new$segments_inferred)))
  cat(sprintf("  -> speedup: %.1fx\n\n", tm_old / tm_new))
}
