# Benchmarks

These scripts compare the package against the original `code/CPLASS.R`
research scripts (cloned separately; not part of this package). They assume
the original repo is checked out at `/home/claude/CPLASS` — adjust the paths
at the top of each script to point at your own checkout.

- `benchmark_correctness_fit.R`: confirms `piecewise_linear_con()` matches
  the original bit-for-bit (to floating point tolerance) on 0/1/2/4-changepoint
  configurations.
- `benchmark_speed_fit.R`: times `piecewise_linear_con()` alone, old vs new,
  across trajectory lengths.
- `benchmark_speed_end_to_end.R`: times full `CPLASS()` runs, old vs new, on
  both real data (`Real_21_Periphery.csv`) and synthetic trajectories.
