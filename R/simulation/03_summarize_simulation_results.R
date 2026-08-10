# 03_summarize_simulation_results.R
# Create compact summary tables corresponding to the main-text simulation tables.
#
# Run this script after:
#   R/simulation/01_study1_exact_equivalence.R
#   R/simulation/02_study2_falsification_diagnostic.R

results_dir <- file.path("results", "simulation")

study1 <- read.csv(
  file.path(results_dir, "study1_condition_results.csv"),
  stringsAsFactors = FALSE
)

study2 <- read.csv(
  file.path(results_dir, "study2_condition_results.csv"),
  stringsAsFactors = FALSE
)


# -------------------------------------------------------------------------
# Study 1 summary:
# For each tau x N combination:
#   - maximum absolute bias across distribution x rho
#   - ranges of SE/SD, coverage, and rejection
# -------------------------------------------------------------------------

study1_groups <- split(
  study1,
  interaction(study1$tau, study1$N, drop = TRUE)
)

study1_summary <- do.call(
  rbind,
  lapply(study1_groups, function(x) {
    data.frame(
      tau = x$tau[1],
      N = x$N[1],
      Max_Abs_Bias = max(abs(x$Bias)),
      SE_SD_Min = min(x$SE_SD),
      SE_SD_Max = max(x$SE_SD),
      Coverage_Min = min(x$Coverage),
      Coverage_Max = max(x$Coverage),
      Rejection_Min = min(x$Rejection),
      Rejection_Max = max(x$Rejection)
    )
  })
)

study1_summary <- study1_summary[
  order(-study1_summary$tau, study1_summary$N),
]
row.names(study1_summary) <- NULL

write.csv(
  study1_summary,
  file.path(results_dir, "study1_main_table_summary.csv"),
  row.names = FALSE
)


# -------------------------------------------------------------------------
# Study 2 summary:
# For each |alpha| x N combination:
#   - maximum absolute bias across rho and signs of alpha
#   - ranges of SE/SD, coverage, and rejection
# -------------------------------------------------------------------------

study2$abs_alpha <- abs(study2$alpha)

study2_groups <- split(
  study2,
  interaction(study2$abs_alpha, study2$N, drop = TRUE)
)

study2_summary <- do.call(
  rbind,
  lapply(study2_groups, function(x) {
    data.frame(
      abs_alpha = x$abs_alpha[1],
      N = x$N[1],
      Max_Abs_Bias = max(abs(x$Bias)),
      SE_SD_Min = min(x$SE_SD),
      SE_SD_Max = max(x$SE_SD),
      Coverage_Min = min(x$Coverage),
      Coverage_Max = max(x$Coverage),
      Rejection_Min = min(x$Rejection),
      Rejection_Max = max(x$Rejection)
    )
  })
)

study2_summary <- study2_summary[
  order(study2_summary$abs_alpha, study2_summary$N),
]
row.names(study2_summary) <- NULL

write.csv(
  study2_summary,
  file.path(results_dir, "study2_main_table_summary.csv"),
  row.names = FALSE
)

cat("\nStudy 1 main-table summary:\n")
print(study1_summary, row.names = FALSE)

cat("\nStudy 2 main-table summary:\n")
print(study2_summary, row.names = FALSE)

cat("\nSaved summary tables in results/simulation/.\n")
