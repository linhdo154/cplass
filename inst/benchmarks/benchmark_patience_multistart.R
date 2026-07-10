suppressMessages(library(cplass))
suppressMessages(library(readr))
suppressMessages(library(dplyr))

data <- suppressMessages(read_csv("/home/claude/CPLASS/data/Real_21_Periphery.csv", show_col_types = FALSE))

for (pid in c(8, 34)) {
  p <- data %>% filter(index_path == pid)
  cat(sprintf("=== path %d (n=%d obs) ===\n", pid, nrow(p)))

  # "gold standard": 5 full-length runs (no early stopping), take the best
  set.seed(1000 + pid)
  tm_gold <- system.time(
    gold <- CPLASS_multistart(p$t, p$x, p$y, iter_max = 5000, burn_in = 300,
                                n_starts = 5, patience = NULL, show_progress = FALSE)
  )["elapsed"]

  # fast+safe: 5 chains with auto patience (0.3 * iter_max = 1500), take the best
  set.seed(1000 + pid)
  tm_fast <- system.time(
    fast <- CPLASS_multistart(p$t, p$x, p$y, iter_max = 5000, burn_in = 300,
                                n_starts = 5, patience = "auto", show_progress = FALSE)
  )["elapsed"]

  cat(sprintf("gold (no patience, 5 starts): best_score=%.4f | n_seg=%d | time=%.2fs\n",
              gold$best$best_score, nrow(gold$best$segments_inferred), tm_gold))
  cat(sprintf("fast (patience=%d, 5 starts): best_score=%.4f | n_seg=%d | time=%.2fs\n",
              fast$patience_used, fast$best$best_score, nrow(fast$best$segments_inferred), tm_fast))
  cat(sprintf("score gap: %.4f (%.3f%%) | time saved: %.1f%%\n\n",
              gold$best$best_score - fast$best$best_score,
              100 * (gold$best$best_score - fast$best$best_score) / abs(gold$best$best_score),
              100 * (1 - tm_fast / tm_gold)))
}
