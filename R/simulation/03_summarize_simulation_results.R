# 03_summarize_simulation_results.R
# Summarize the Study 1 and Study 2 results in the CURRENT manuscript format.
#
# Run this script after:
#   R/simulation/01_study1_exact_equivalence.R
#   R/simulation/02_study2_falsification_diagnostic.R
#
# No random numbers are generated in this script.

source(file.path("R", "simulation", "00_helpers.R"))

study1_file <- file.path(SETTINGS$output_dir, "study1_condition_results.csv")
study2_file <- file.path(SETTINGS$output_dir, "study2_condition_results.csv")

if (!file.exists(study1_file)) {
  stop("Study 1 results not found. Run 01_study1_exact_equivalence.R first.")
}
if (!file.exists(study2_file)) {
  stop("Study 2 results not found. Run 02_study2_falsification_diagnostic.R first.")
}

study1 <- read.csv(study1_file, stringsAsFactors = FALSE, check.names = FALSE)
study2 <- read.csv(study2_file, stringsAsFactors = FALSE, check.names = FALSE)

# -----------------------------------------------------------------------------
# Study 1: current main-text table, summarized by tau x N
# -----------------------------------------------------------------------------

study1_main_rows <- list()
k <- 1L

for (tau_value in TAU_LEVELS) {
  for (n_value in N_LEVELS) {
    x <- study1[
      abs(study1$tau - tau_value) < 1e-12 & study1$N == n_value,
    ]

    study1_main_rows[[k]] <- data.frame(
      tau = tau_value,
      N = n_value,
      Max_Abs_Bias = max(abs(x$mc_bias)),
      SE_SD_Min = min(x$se_sd),
      SE_SD_Max = max(x$se_sd),
      Coverage_Min = min(x$coverage),
      Coverage_Max = max(x$coverage),
      Rejection_Min = min(x$reject_beta),
      Rejection_Max = max(x$reject_beta),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

study1_main <- do.call(rbind, study1_main_rows)
row.names(study1_main) <- NULL

write.csv(
  study1_main,
  file.path(SETTINGS$output_dir, "study1_main_table_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Study 2: current main-text table, summarized by |alpha| x N
# -----------------------------------------------------------------------------

alpha_magnitudes <- c(0, 0.025, 0.05, 0.075, 0.10)
study2_main_rows <- list()
k <- 1L

for (alpha_abs in alpha_magnitudes) {
  for (n_value in N_LEVELS) {
    x <- study2[
      abs(abs(study2$alpha) - alpha_abs) < 1e-12 & study2$N == n_value,
    ]

    study2_main_rows[[k]] <- data.frame(
      abs_alpha = alpha_abs,
      N = n_value,
      Max_Abs_Bias = max(abs(x$mc_bias_alpha)),
      SE_SD_Min = min(x$se_sd_alpha),
      SE_SD_Max = max(x$se_sd_alpha),
      Coverage_Min = min(x$coverage_alpha),
      Coverage_Max = max(x$coverage_alpha),
      Rejection_Min = min(x$reject_alpha),
      Rejection_Max = max(x$reject_alpha),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

study2_main <- do.call(rbind, study2_main_rows)
row.names(study2_main) <- NULL

write.csv(
  study2_main,
  file.path(SETTINGS$output_dir, "study2_main_table_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Appendix-facing compact files in the same order as the manuscript tables
# -----------------------------------------------------------------------------

study1_appendix_compact <- study1[, c(
  "N", "tau", "rho", "distribution", "mc_bias", "se_sd", "coverage", "reject_beta"
)]

study1_appendix_compact <- study1_appendix_compact[order(
  study1_appendix_compact$N,
  match(study1_appendix_compact$tau, TAU_LEVELS),
  study1_appendix_compact$rho,
  match(study1_appendix_compact$distribution, DIST_LEVELS)
), ]

names(study1_appendix_compact) <- c(
  "N", "tau", "rho", "Distribution", "Bias", "SE_SD", "Coverage", "Rejection"
)
row.names(study1_appendix_compact) <- NULL

write.csv(
  study1_appendix_compact,
  file.path(SETTINGS$output_dir, "study1_appendix_manuscript.csv"),
  row.names = FALSE
)

study2_appendix_compact <- study2[, c(
  "N", "alpha", "rho", "mc_bias_alpha", "se_sd_alpha", "coverage_alpha", "reject_alpha"
)]

study2_appendix_compact <- study2_appendix_compact[order(
  study2_appendix_compact$N,
  match(study2_appendix_compact$alpha, ALPHA_LEVELS),
  study2_appendix_compact$rho
), ]

names(study2_appendix_compact) <- c(
  "N", "alpha", "rho", "Bias", "SE_SD", "Coverage", "Rejection"
)
row.names(study2_appendix_compact) <- NULL

write.csv(
  study2_appendix_compact,
  file.path(SETTINGS$output_dir, "study2_appendix_manuscript.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Automatic check against the values currently printed in the manuscript
# -----------------------------------------------------------------------------

expected_study1 <- data.frame(
  tau = c(0, 0, 0, -0.05, -0.05, -0.05, -0.10, -0.10, -0.10),
  N = c(500, 1000, 2000, 500, 1000, 2000, 500, 1000, 2000),
  Max_Abs_Bias = c(.0016, .0005, .0006, .0015, .0010, .0005, .0015, .0009, .0006),
  SE_SD_Min = c(.991, .991, .982, .982, .986, .966, .992, .986, .977),
  SE_SD_Max = c(1.017, 1.018, 1.010, 1.019, 1.017, 1.007, 1.013, 1.007, 1.016),
  Coverage_Min = c(.944, .946, .944, .946, .943, .941, .944, .943, .941),
  Coverage_Max = c(.953, .953, .954, .955, .956, .951, .955, .954, .955),
  Rejection_Min = c(.047, .047, .046, .141, .229, .421, .421, .701, .941),
  Rejection_Max = c(.056, .054, .056, .210, .367, .639, .626, .893, .996)
)

expected_study2 <- data.frame(
  abs_alpha = rep(c(0, .025, .05, .075, .10), each = 3),
  N = rep(c(500, 1000, 2000), times = 5),
  Max_Abs_Bias = c(
    .0005, .0005, .0003,
    .0006, .0007, .0003,
    .0005, .0007, .0004,
    .0011, .0002, .0005,
    .0005, .0004, .0004
  ),
  SE_SD_Min = c(
    .971, .997, .977,
    .985, .981, .984,
    .974, .991, .993,
    .980, .996, .989,
    .980, .982, .982
  ),
  SE_SD_Max = c(
    1.002, 1.012, 1.021,
    1.020, 1.023, 1.021,
    1.010, 1.013, 1.005,
    1.016, 1.004, 1.004,
    1.001, 1.009, 1.007
  ),
  Coverage_Min = c(
    .947, .950, .946,
    .945, .943, .945,
    .944, .947, .947,
    .941, .946, .946,
    .945, .946, .943
  ),
  Coverage_Max = c(
    .950, .951, .957,
    .953, .954, .954,
    .954, .954, .954,
    .957, .952, .951,
    .952, .952, .953
  ),
  Rejection_Min = c(
    .050, .049, .043,
    .094, .143, .235,
    .226, .412, .686,
    .442, .748, .955,
    .683, .937, .998
  ),
  Rejection_Max = c(
    .053, .050, .054,
    .124, .202, .355,
    .353, .607, .883,
    .659, .919, .998,
    .885, .995, 1.000
  )
)

format_main <- function(x) {
  data.frame(
    key1 = sprintf("%.3f", x[[1]]),
    N = as.integer(x$N),
    Max_Abs_Bias = sprintf("%.4f", x$Max_Abs_Bias),
    SE_SD_Min = sprintf("%.3f", x$SE_SD_Min),
    SE_SD_Max = sprintf("%.3f", x$SE_SD_Max),
    Coverage_Min = sprintf("%.3f", x$Coverage_Min),
    Coverage_Max = sprintf("%.3f", x$Coverage_Max),
    Rejection_Min = sprintf("%.3f", x$Rejection_Min),
    Rejection_Max = sprintf("%.3f", x$Rejection_Max),
    stringsAsFactors = FALSE
  )
}

study1_pass <- identical(format_main(study1_main), format_main(expected_study1))
study2_pass <- identical(format_main(study2_main), format_main(expected_study2))

study1_max_abs_bias <- max(abs(study1$mc_bias))
study1_se_sd_range <- range(study1$se_sd)
study1_coverage_range <- range(study1$coverage)
study1_type1_range <- range(study1$reject_beta[abs(study1$tau) < 1e-12])

study2_max_abs_bias <- max(abs(study2$mc_bias_alpha))
study2_se_sd_range <- range(study2$se_sd_alpha)
study2_coverage_range <- range(study2$coverage_alpha)
study2_type1_range <- range(study2$reject_alpha[abs(study2$alpha) < 1e-12])

check_lines <- c(
  paste0("Study 1 main-text table: ", if (study1_pass) "PASS" else "FAIL"),
  paste0("Study 2 main-text table: ", if (study2_pass) "PASS" else "FAIL"),
  "",
  paste0("Study 1 max absolute bias: ", sprintf("%.4f", study1_max_abs_bias)),
  paste0(
    "Study 1 SE/SD range: ",
    sprintf("%.3f", study1_se_sd_range[1]), " to ",
    sprintf("%.3f", study1_se_sd_range[2])
  ),
  paste0(
    "Study 1 coverage range: ",
    sprintf("%.3f", study1_coverage_range[1]), " to ",
    sprintf("%.3f", study1_coverage_range[2])
  ),
  paste0(
    "Study 1 Type I range: ",
    sprintf("%.3f", study1_type1_range[1]), " to ",
    sprintf("%.3f", study1_type1_range[2])
  ),
  "",
  paste0("Study 2 max absolute bias: ", sprintf("%.4f", study2_max_abs_bias)),
  paste0(
    "Study 2 SE/SD range: ",
    sprintf("%.3f", study2_se_sd_range[1]), " to ",
    sprintf("%.3f", study2_se_sd_range[2])
  ),
  paste0(
    "Study 2 coverage range: ",
    sprintf("%.3f", study2_coverage_range[1]), " to ",
    sprintf("%.3f", study2_coverage_range[2])
  ),
  paste0(
    "Study 2 Type I range: ",
    sprintf("%.3f", study2_type1_range[1]), " to ",
    sprintf("%.3f", study2_type1_range[2])
  )
)

writeLines(
  check_lines,
  file.path(SETTINGS$output_dir, "manuscript_reproduction_check.txt")
)

cat("\n", paste(check_lines, collapse = "\n"), "\n", sep = "")
cat("\nSummary and appendix-facing files saved to: ", SETTINGS$output_dir, "\n", sep = "")
