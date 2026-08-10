# 02_study2_falsification_diagnostic.R
# Study 2: Estimation and Inference for the Falsification Diagnostic
#
# Reproduces Study 2 results reported in the manuscript.
# Run from the repository root.

source(file.path("R", "simulation", "00_helpers.R"))

cat("Preparing Study 2 conditions...\n")

# -----------------------------------------------------------------------------
# Calibrate c(alpha)
# -----------------------------------------------------------------------------

shift_table <- data.frame(
  alpha = ALPHA_LEVELS,
  shift_c = vapply(ALPHA_LEVELS, find_anchor_shift, numeric(1)),
  stringsAsFactors = FALSE
)

shift_table$alpha_check <- vapply(
  shift_table$shift_c,
  alpha_from_shift,
  numeric(1)
)

if (max(abs(shift_table$alpha_check - shift_table$alpha)) > 1e-8) {
  stop("The anchor-shift calibration did not recover the target alpha values.")
}

study2_calibration_out <- data.frame(
  alpha = shift_table$alpha,
  c_alpha = shift_table$shift_c,
  stringsAsFactors = FALSE
)

write.csv(
  study2_calibration_out,
  file.path(SETTINGS$output_dir, "study2_calibrated_c_alpha.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Build the 81 simulation conditions
# -----------------------------------------------------------------------------

study2_grid <- expand.grid(
  N = N_LEVELS,
  rho = RHO_LEVELS,
  mu1 = MU1_DIAGNOSTIC,
  alpha = ALPHA_LEVELS,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

study2_grid <- merge(
  study2_grid,
  shift_table[, c("alpha", "shift_c")],
  by = "alpha",
  sort = FALSE
)

study2_grid <- study2_grid[order(
  study2_grid$N,
  study2_grid$rho,
  match(study2_grid$alpha, ALPHA_LEVELS)
), ]

row.names(study2_grid) <- NULL
study2_grid$condition_id <- seq_len(nrow(study2_grid))

study2_grid$condition_seed <-
  SETTINGS$seed + 500000L + study2_grid$condition_id * 1013L

if (nrow(study2_grid) != 81L) {
  stop("Study 2 should contain exactly 81 conditions; found ", nrow(study2_grid), ".")
}

# -----------------------------------------------------------------------------
# Run Study 2
# -----------------------------------------------------------------------------

cat(
  "Running Study 2:", nrow(study2_grid), "conditions x",
  SETTINGS$n_rep, "replications.\n"
)

study2_results <- do.call(
  rbind,
  parallel_lapply_conditions(
    split(study2_grid, seq_len(nrow(study2_grid))),
    simulate_study2_condition
  )
)

row.names(study2_results) <- NULL
study2_results <- study2_results[order(
  study2_results$N,
  study2_results$rho,
  match(study2_results$alpha, ALPHA_LEVELS)
), ]

write.csv(
  study2_results,
  file.path(SETTINGS$output_dir, "study2_condition_results.csv"),
  row.names = FALSE
)

# Appendix-ready complete results.
study2_appendix_full <- study2_results[, c(
  "alpha", "N", "rho", "mu1", "shift_c", "replications",
  "mc_mean_alpha", "mc_bias_alpha", "mean_se_alpha",
  "empirical_sd_alpha", "se_sd_alpha", "coverage_alpha", "reject_alpha"
)]

write.csv(
  study2_appendix_full,
  file.path(SETTINGS$output_dir, "study2_appendix_full.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(SETTINGS$output_dir, "sessionInfo_study2.txt")
)

cat("\nStudy 2 complete.\n")
cat("Saved to: ", SETTINGS$output_dir, "\n", sep = "")
