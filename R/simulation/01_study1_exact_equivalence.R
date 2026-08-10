# 01_study1_exact_equivalence.R
# Study 1: Finite-Sample Inference Under Exact Identification
#
# Reproduces Study 1 results reported in the manuscript.
# Run from the repository root.

source(file.path("R", "simulation", "00_helpers.R"))

cat("Preparing Study 1 conditions...\n")

# -----------------------------------------------------------------------------
# Calibrate delta_T
# -----------------------------------------------------------------------------

study1_delta_table <- expand.grid(
  distribution = DIST_LEVELS,
  tau = TAU_LEVELS,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

study1_delta_table$delta_t <- mapply(
  FUN = find_delta_t,
  target_tau = study1_delta_table$tau,
  distribution = study1_delta_table$distribution
)

study1_delta_table$tau_check <- mapply(
  FUN = function(distribution, delta_t) {
    integrate_safe(function(theta) {
      (p_test_focal(theta, delta_t) - p_test_reference(theta)) *
        f1_density_study1(theta, distribution)
    })
  },
  distribution = study1_delta_table$distribution,
  delta_t = study1_delta_table$delta_t
)

if (max(abs(study1_delta_table$tau_check - study1_delta_table$tau)) > 1e-8) {
  stop("The Study 1 delta_t calibration did not recover the target tau values.")
}

# A compact calibration file matching the manuscript notation.
study1_calibration_out <- data.frame(
  Distribution = study1_delta_table$distribution,
  tau = study1_delta_table$tau,
  delta_T = study1_delta_table$delta_t,
  stringsAsFactors = FALSE
)

write.csv(
  study1_calibration_out,
  file.path(SETTINGS$output_dir, "study1_calibrated_delta_T.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Build the 81 simulation conditions
# -----------------------------------------------------------------------------

study1_grid <- expand.grid(
  N = N_LEVELS,
  rho = RHO_LEVELS,
  distribution = DIST_LEVELS,
  tau = TAU_LEVELS,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

study1_grid <- merge(
  study1_grid,
  study1_delta_table[, c("distribution", "tau", "delta_t")],
  by = c("distribution", "tau"),
  sort = FALSE
)

study1_grid <- study1_grid[order(
  study1_grid$N,
  study1_grid$rho,
  match(study1_grid$distribution, DIST_LEVELS),
  match(study1_grid$tau, TAU_LEVELS)
), ]

row.names(study1_grid) <- NULL
study1_grid$condition_id <- seq_len(nrow(study1_grid))

study1_grid$condition_seed <-
  SETTINGS$seed + study1_grid$condition_id * 1013L

if (nrow(study1_grid) != 81L) {
  stop("Study 1 should contain exactly 81 conditions; found ", nrow(study1_grid), ".")
}

# -----------------------------------------------------------------------------
# Run Study 1
# -----------------------------------------------------------------------------

cat(
  "Running Study 1:", nrow(study1_grid), "conditions x",
  SETTINGS$n_rep, "replications.\n"
)

study1_results <- do.call(
  rbind,
  parallel_lapply_conditions(
    split(study1_grid, seq_len(nrow(study1_grid))),
    simulate_study1_condition
  )
)

row.names(study1_results) <- NULL
study1_results <- study1_results[order(
  study1_results$N,
  study1_results$rho,
  match(study1_results$distribution, DIST_LEVELS),
  match(study1_results$tau, TAU_LEVELS)
), ]

write.csv(
  study1_results,
  file.path(SETTINGS$output_dir, "study1_condition_results.csv"),
  row.names = FALSE
)

# Appendix-ready complete results.
study1_appendix_full <- study1_results[, c(
  "N", "rho", "distribution", "tau", "delta_t", "replications",
  "mc_mean_did", "mc_bias", "rmse", "mean_se", "empirical_sd",
  "se_sd", "coverage", "reject_beta", "mean_alpha",
  "mean_se_alpha", "empirical_sd_alpha", "reject_alpha", "coverage_alpha"
)]

names(study1_appendix_full)[names(study1_appendix_full) == "reject_beta"] <-
  "rejection_rate_beta"

write.csv(
  study1_appendix_full,
  file.path(SETTINGS$output_dir, "study1_appendix_full.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(SETTINGS$output_dir, "sessionInfo_study1.txt")
)

cat("\nStudy 1 complete.\n")
cat("Saved to: ", SETTINGS$output_dir, "\n", sep = "")
