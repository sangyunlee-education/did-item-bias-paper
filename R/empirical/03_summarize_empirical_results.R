# 03_summarize_empirical_results.R
# Format empirical results and check them against the manuscript.
#
# Run this script after:
#   R/empirical/01_conventional_dif_screening.R
#   R/empirical/02_pairwise_falsification_did.R

results_dir <- file.path("results", "empirical")

DATA_FILE <- file.path(
  "data", "derived", "piaac_korea_locator_analysis.csv"
)

SCREENING_FILE <- file.path(
  results_dir,
  "conventional_dif_screening.csv"
)

PAIR_FILE <- file.path(
  results_dir,
  "pairwise_falsification_did.csv"
)

required_files <- c(
  DATA_FILE,
  SCREENING_FILE,
  PAIR_FILE
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required output files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

dat <- utils::read.csv(
  DATA_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

screening <- utils::read.csv(
  SCREENING_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pairs <- utils::read.csv(
  PAIR_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# Formatting helpers
# ------------------------------------------------------------------------------

format_p <- function(p) {
  ifelse(
    p < 0.001,
    "<.001",
    sub("^0", "", sprintf("%.3f", p))
  )
}

format_num3 <- function(x) {
  out <- sprintf("%.3f", x)
  out <- sub("^0", "", out)
  out <- sub("^-0", "-", out)
  out
}

format_est_se <- function(estimate, se) {
  paste0(
    format_num3(estimate),
    " (",
    format_num3(se),
    ")"
  )
}

format_interval <- function(lower, upper) {
  paste0(
    "[",
    format_num3(lower),
    ", ",
    format_num3(upper),
    "]"
  )
}


# ------------------------------------------------------------------------------
# Formatted Appendix screening table
# ------------------------------------------------------------------------------

screening_table <- data.frame(
  Item = screening$Item,
  `Mantel--Haenszel` = format_p(screening$Mantel_Haenszel),
  `Logistic DIF` = format_p(screening$Logistic_DIF),
  SIBTEST = format_p(screening$SIBTEST),
  `2PL IRT-LRT` = format_p(screening$IRT_2PL_LRT),
  `Candidate anchor` = ifelse(
    screening$Candidate_anchor,
    "Yes",
    "No"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

utils::write.csv(
  screening_table,
  file.path(results_dir, "screening_table_formatted.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Formatted main-text Table 3: all 15 test--anchor pairs
# ------------------------------------------------------------------------------

pairwise_table <- data.frame(
  `Test item` = pairs$Test_item,
  `Anchor item` = pairs$Anchor_item,
  `alpha (SE)` = format_est_se(
    pairs$alpha,
    pairs$SE_alpha
  ),
  `p_alpha` = format_p(pairs$p_alpha),
  Diagnostic = pairs$Diagnostic,
  `DID (SE)` = format_est_se(
    pairs$DID,
    pairs$SE_DID
  ),
  `95% CI` = format_interval(
    pairs$CI_lower,
    pairs$CI_upper
  ),
  `p_DID` = format_p(pairs$p_DID),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

utils::write.csv(
  pairwise_table,
  file.path(results_dir, "pairwise_table_formatted.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# Manuscript reproduction checks
# ------------------------------------------------------------------------------

expected_candidate_anchors <- c(
  "C601C06S",
  "C815P001S",
  "C833P001S"
)

expected_screening <- data.frame(
  Item = c(
    "C601C06S",
    "C815P001S",
    "C815P002S",
    "C832P001S",
    "C832P002S",
    "C813P001S",
    "C833P001S",
    "C833P002S"
  ),
  Mantel_Haenszel = c(
    ".159", ".162", ".328", ".365",
    ".016", "<.001", ".668", ".003"
  ),
  Logistic_DIF = c(
    ".113", ".109", ".099", ".236",
    ".037", "<.001", ".747", ".005"
  ),
  SIBTEST = c(
    ".434", ".968", ".002", ".218",
    ".687", ".022", ".335", ".268"
  ),
  IRT_2PL_LRT = c(
    ".083", ".377", ".032", ".005",
    ".561", ".006", ".121", ".176"
  ),
  Candidate_anchor = c(
    "Yes", "Yes", "No", "No",
    "No", "No", "Yes", "No"
  ),
  stringsAsFactors = FALSE
)

observed_screening <- data.frame(
  Item = screening$Item,
  Mantel_Haenszel = format_p(screening$Mantel_Haenszel),
  Logistic_DIF = format_p(screening$Logistic_DIF),
  SIBTEST = format_p(screening$SIBTEST),
  IRT_2PL_LRT = format_p(screening$IRT_2PL_LRT),
  Candidate_anchor = ifelse(
    screening$Candidate_anchor,
    "Yes",
    "No"
  ),
  stringsAsFactors = FALSE
)

expected_pairs <- data.frame(
  Test_item = c(
    "C815P002S", "C815P002S", "C815P002S",
    "C832P001S", "C832P001S", "C832P001S",
    "C832P002S", "C832P002S", "C832P002S",
    "C813P001S", "C813P001S", "C813P001S",
    "C833P002S", "C833P002S", "C833P002S"
  ),
  Anchor_item = rep(
    c("C601C06S", "C815P001S", "C833P001S"),
    5
  ),
  alpha = c(
    "-.088", "-.104", ".071",
    "-.148", "-.164", ".011",
    "-.134", "-.150", ".025",
    "-.147", "-.163", ".013",
    "-.356", "-.372", "-.197"
  ),
  SE_alpha = c(
    ".008", ".007", ".009",
    ".009", ".008", ".010",
    ".009", ".008", ".010",
    ".009", ".009", ".010",
    ".010", ".010", ".011"
  ),
  p_alpha = c(
    "<.001", "<.001", "<.001",
    "<.001", "<.001", ".265",
    "<.001", "<.001", ".010",
    "<.001", "<.001", ".231",
    "<.001", "<.001", "<.001"
  ),
  Diagnostic = c(
    "Rejected", "Rejected", "Rejected",
    "Rejected", "Rejected", "Did not reject",
    "Rejected", "Rejected", "Rejected",
    "Rejected", "Rejected", "Did not reject",
    "Rejected", "Rejected", "Rejected"
  ),
  DID = c(
    ".009", ".008", ".017",
    ".010", ".010", ".019",
    "-.026", "-.027", "-.018",
    "-.033", "-.034", "-.025",
    "-.037", "-.038", "-.029"
  ),
  SE_DID = c(
    ".011", ".010", ".013",
    ".012", ".012", ".014",
    ".012", ".012", ".014",
    ".013", ".013", ".015",
    ".014", ".014", ".016"
  ),
  CI_lower = c(
    "-.013", "-.013", "-.009",
    "-.014", "-.014", "-.009",
    "-.050", "-.050", "-.044",
    "-.058", "-.058", "-.054",
    "-.065", "-.065", "-.059"
  ),
  CI_upper = c(
    ".030", ".028", ".043",
    ".034", ".033", ".046",
    "-.003", "-.004", ".009",
    "-.008", "-.009", ".004",
    "-.009", "-.010", ".002"
  ),
  p_DID = c(
    ".431", ".446", ".197",
    ".397", ".417", ".184",
    ".030", ".020", ".187",
    ".009", ".007", ".097",
    ".009", ".007", ".065"
  ),
  stringsAsFactors = FALSE
)

observed_pairs <- data.frame(
  Test_item = pairs$Test_item,
  Anchor_item = pairs$Anchor_item,
  alpha = format_num3(pairs$alpha),
  SE_alpha = format_num3(pairs$SE_alpha),
  p_alpha = format_p(pairs$p_alpha),
  Diagnostic = pairs$Diagnostic,
  DID = format_num3(pairs$DID),
  SE_DID = format_num3(pairs$SE_DID),
  CI_lower = format_num3(pairs$CI_lower),
  CI_upper = format_num3(pairs$CI_upper),
  p_DID = format_p(pairs$p_DID),
  stringsAsFactors = FALSE
)

sample_pass <- (
  nrow(dat) == 5780L &&
  sum(dat$G == 0L) == 2715L &&
  sum(dat$G == 1L) == 3065L
)

candidate_anchor_pass <- identical(
  screening$Item[screening$Candidate_anchor],
  expected_candidate_anchors
)

screening_pass <- identical(
  observed_screening,
  expected_screening
)

pairwise_pass <- identical(
  observed_pairs,
  expected_pairs
)

rejected <- pairs[pairs$Diagnostic == "Rejected", , drop = FALSE]
nonrejected <- pairs[pairs$Diagnostic == "Did not reject", , drop = FALSE]

narrative_pass <- (
  nrow(pairs) == 15L &&
  nrow(rejected) == 13L &&
  nrow(nonrejected) == 2L &&
  sum(rejected$p_DID < 0.05) == 6L
)

checks <- data.frame(
  Check = c(
    "Analytic sample",
    "Candidate anchors",
    "Appendix DIF screening table",
    "Main-text Table 3",
    "Narrative counts"
  ),
  Status = c(
    ifelse(sample_pass, "PASS", "FAIL"),
    ifelse(candidate_anchor_pass, "PASS", "FAIL"),
    ifelse(screening_pass, "PASS", "FAIL"),
    ifelse(pairwise_pass, "PASS", "FAIL"),
    ifelse(narrative_pass, "PASS", "FAIL")
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  checks,
  file.path(
    results_dir,
    "manuscript_reproduction_check.csv"
  ),
  row.names = FALSE
)

cat("\nManuscript reproduction checks\n")
cat("------------------------------\n")
print(checks, row.names = FALSE)

overall_pass <- all(checks$Status == "PASS")

cat("\nOVERALL: ", ifelse(overall_pass, "PASS", "FAIL"), "\n", sep = "")

if (!screening_pass) {
  screening_compare <- cbind(
    expected_screening,
    observed_screening[
      setdiff(names(observed_screening), "Item")
    ]
  )

  utils::write.csv(
    screening_compare,
    file.path(
      results_dir,
      "screening_manuscript_comparison.csv"
    ),
    row.names = FALSE
  )

  cat(
    "\nScreening mismatch detected. See ",
    "results/empirical/screening_manuscript_comparison.csv\n",
    sep = ""
  )
}

if (!pairwise_pass) {
  pairwise_compare <- cbind(
    expected_pairs,
    observed_pairs[
      setdiff(names(observed_pairs), c("Test_item", "Anchor_item"))
    ]
  )

  utils::write.csv(
    pairwise_compare,
    file.path(
      results_dir,
      "pairwise_manuscript_comparison.csv"
    ),
    row.names = FALSE
  )

  cat(
    "\nTable 3 mismatch detected. See ",
    "results/empirical/pairwise_manuscript_comparison.csv\n",
    sep = ""
  )
}

cat("\nFormatted tables saved in results/empirical/.\n")
