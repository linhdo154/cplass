# Package index

## Main algorithm

Run CPLASS and manage the Metropolis-Hastings search.

- [`CPLASS()`](https://linhdo154.github.io/cplass/reference/CPLASS-function.md)
  : Run CPLASS on a single 2D trajectory
- [`CPLASS_multistart()`](https://linhdo154.github.io/cplass/reference/CPLASS_multistart.md)
  : Run CPLASS from multiple independent random starts and keep the best
- [`MHsearch()`](https://linhdo154.github.io/cplass/reference/MHsearch.md)
  : Metropolis-Hastings search over the changepoint space
- [`FinalMH()`](https://linhdo154.github.io/cplass/reference/FinalMH.md)
  : Collect a Metropolis-Hastings trace into flat vectors/lists

## Convergence diagnostics

Check whether the stochastic search has stabilized.

- [`check_convergence()`](https://linhdo154.github.io/cplass/reference/check_convergence.md)
  : Check whether a single CPLASS chain has stopped improving
- [`plot_convergence()`](https://linhdo154.github.io/cplass/reference/plot_convergence.md)
  : Plot a CPLASS chain's convergence trace

## Piecewise-linear fitting

Continuous piecewise-linear MLE given a set of changepoints (Section 2.2
of the paper).

- [`piecewise_linear_con()`](https://linhdo154.github.io/cplass/reference/piecewise_linear_con.md)
  : Continuous piecewise-linear fit given a set of changepoints
- [`build_j_inferred()`](https://linhdo154.github.io/cplass/reference/build_j_inferred.md)
  : Broadcast per-segment states back onto the original observation grid

## Criterion function and penalties

The penalized-likelihood criterion Phi(r) (Section 2.3).

- [`CS()`](https://linhdo154.github.io/cplass/reference/CS.md) :
  Criterion function for a changepoint configuration
- [`loglikelihood()`](https://linhdo154.github.io/cplass/reference/loglikelihood.md)
  : Gaussian log-likelihood from residual sum of squares
- [`sSIC()`](https://linhdo154.github.io/cplass/reference/sSIC.md) :
  Strengthened Schwarz Information Criterion penalty
- [`AICc()`](https://linhdo154.github.io/cplass/reference/AICc.md) :
  Corrected AIC penalty
- [`cp_to_r()`](https://linhdo154.github.io/cplass/reference/cp_to_r.md)
  : Convert a changepoint index vector to a 0/1 indicator vector
- [`r_to_cp()`](https://linhdo154.github.io/cplass/reference/r_to_cp.md)
  : Convert a 0/1 changepoint indicator vector to changepoint indices

## MCMC proposal mechanisms

The four changepoint-proposal types described in Section 2.4.

- [`log_q_new_pmf()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`q_new()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`log_q_bd_pmf()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`q_bd()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`q_shift()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`q_as_pmf()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`q_as()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`q_ds_pmf()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`q_ds()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`q_bd2()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`log_q_bd2_pmf()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`proposal_function()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  [`pproposal()`](https://linhdo154.github.io/cplass/reference/proposals.md)
  : Proposal functions for the CPLASS Metropolis-Hastings sampler

## Cumulative Speed Allocation (CSA)

Duration-weighted speed summaries (Cook et al. 2025).

- [`compute_csa()`](https://linhdo154.github.io/cplass/reference/compute_csa.md)
  : Cumulative Speed Allocation (CSA)
- [`csa_theoretical()`](https://linhdo154.github.io/cplass/reference/csa_theoretical.md)
  : Theoretical CSA under a two-state (stationary/motile) kinetic model
- [`plot_csa()`](https://linhdo154.github.io/cplass/reference/plot_csa.md)
  : Plot Cumulative Speed Allocation curves
- [`summarize_segments()`](https://linhdo154.github.io/cplass/reference/summarize_segments.md)
  : Stack segment summaries from a list of true (ground-truth) paths
- [`summarize_segments_inferred()`](https://linhdo154.github.io/cplass/reference/summarize_segments_inferred.md)
  : Stack segment summaries from a list of CPLASS-inferred paths
- [`summarize_segments_cutoff()`](https://linhdo154.github.io/cplass/reference/summarize_segments_cutoff.md)
  : Stack segment summaries from inferred paths (cutoff-labeled states)
- [`summarize_segments_test()`](https://linhdo154.github.io/cplass/reference/summarize_segments_test.md)
  : Stack segment summaries using a threshold-based test-state column

## Visualization

- [`plot_path_inferred()`](https://linhdo154.github.io/cplass/reference/plot_path_inferred.md)
  : Four-panel dashboard for a single CPLASS-segmented trajectory
- [`plot_path_inferred_second_version()`](https://linhdo154.github.io/cplass/reference/plot_path_inferred_second_version.md)
  : Alternate 3-row layout for a single segmented trajectory
- [`plot_path_inferred_xy()`](https://linhdo154.github.io/cplass/reference/plot_path_inferred_xy.md)
  : Just the x-y panel of a segmented trajectory
- [`plot_path_actual_and_inferred()`](https://linhdo154.github.io/cplass/reference/plot_path_actual_and_inferred.md)
  : Side-by-side comparison of true vs. inferred segmentation
- [`plot_msd()`](https://linhdo154.github.io/cplass/reference/plot_msd.md)
  : Plot pathwise and ensemble MSD curves for a group, split by activity
- [`get_plot_size()`](https://linhdo154.github.io/cplass/reference/get_plot_size.md)
  : Compute shared plot bounds across a list of ground-truth paths
- [`get_plot_size_inferred()`](https://linhdo154.github.io/cplass/reference/get_plot_size_inferred.md)
  : Compute shared plot bounds across a list of inferred paths
- [`get_one_plot_size()`](https://linhdo154.github.io/cplass/reference/get_one_plot_size.md)
  : Compute plot bounds for a single inferred path

## Mean squared displacement

- [`pathwise_MSD()`](https://linhdo154.github.io/cplass/reference/pathwise_MSD.md)
  : Pathwise mean squared displacement (MSD)
- [`ensemble_MSD()`](https://linhdo154.github.io/cplass/reference/ensemble_MSD.md)
  : Ensemble (group-averaged) MSD

## Utilities

- [`fill_missing_data()`](https://linhdo154.github.io/cplass/reference/fill_missing_data.md)
  : Fill gaps in an unevenly-sampled 2D trajectory
- [`duration_of_difference()`](https://linhdo154.github.io/cplass/reference/duration_of_difference.md)
  : Total duration over which true and inferred segmentations disagree
  on state
- [`infer_states_speed_cutoff()`](https://linhdo154.github.io/cplass/reference/infer_states_speed_cutoff.md)
  : Classify inferred segments as motile/stationary by a speed cutoff
- [`cplass`](https://linhdo154.github.io/cplass/reference/cplass-package.md)
  [`cplass-package`](https://linhdo154.github.io/cplass/reference/cplass-package.md)
  : cplass: Continuous Piecewise-Linear Approximation via Stochastic
  Search
