# 03_summarize_empirical_results.R
# Create compact manuscript-facing summaries for the empirical illustration.
#
# Run this script after:
#   R/empirical/01_conventional_dif_screening.R
#   R/empirical/02_pairwise_falsification_did.R

results_dir <- file.path("results", "empirical")

screening <- read.csv(
  file.path(results_dir, "conventional_dif_screening.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pairs <- read.csv(
  file.path(results_dir, "pairwise_falsification_did.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# -------------------------------------------------------------------------
# Formatting helpers
# -------------------------------------------------------------------------

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


# -------------------------------------------------------------------------
# Appendix screening table
# -------------------------------------------------------------------------

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

write.csv(
  screening_table,
  file.path(results_dir, "screening_table_formatted.csv"),
  row.names = FALSE
)


# -------------------------------------------------------------------------
# Main-text pairwise table
# -------------------------------------------------------------------------

pairwise_table <- data.frame(
  `Test item` = pairs$Test_item,
  `Anchor item` = pairs$Anchor_item,
  `alpha (SE)` = paste0(
    format_num3(pairs$alpha),
    " (",
    format_num3(pairs$SE_alpha),
    ")"
  ),
  `p_alpha` = format_p(pairs$p_alpha),
  Diagnostic = pairs$Diagnostic,
  `DID (SE)` = paste0(
    format_num3(pairs$DID),
    " (",
    format_num3(pairs$SE_DID),
    ")"
  ),
  `95% CI` = paste0(
    "[",
    format_num3(pairs$CI_lower),
    ", ",
    format_num3(pairs$CI_upper),
    "]"
  ),
  `p_DID` = format_p(pairs$p_DID),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  pairwise_table,
  file.path(results_dir, "pairwise_table_formatted.csv"),
  row.names = FALSE
)


# -------------------------------------------------------------------------
# Key counts for checking the manuscript narrative
# -------------------------------------------------------------------------

candidate_anchors <- screening$Item[
  screening$Candidate_anchor
]

nonrejected <- pairs[
  pairs$Diagnostic == "Did not reject",
]

rejected <- pairs[
  pairs$Diagnostic == "Rejected",
]

significant_rejected_did <- rejected[
  rejected$p_DID < 0.05,
]

cat("\nEmpirical illustration summary\n")
cat("------------------------------\n")
cat(
  "Candidate anchors (", length(candidate_anchors), "): ",
  paste(candidate_anchors, collapse = ", "),
  "\n",
  sep = ""
)
cat("Total test--anchor pairs: ", nrow(pairs), "\n", sep = "")
cat("Diagnostic rejected:      ", nrow(rejected), "\n", sep = "")
cat("Diagnostic did not reject:", nrow(nonrejected), "\n", sep = "")
cat(
  "Significant DID among rejected pairs: ",
  nrow(significant_rejected_did),
  "\n",
  sep = ""
)

cat("\nNonrejected pairs:\n")
print(
  nonrejected[
    c(
      "Test_item",
      "Anchor_item",
      "alpha",
      "SE_alpha",
      "p_alpha",
      "DID",
      "SE_DID",
      "CI_lower",
      "CI_upper",
      "p_DID"
    )
  ],
  row.names = FALSE
)

cat("\nSaved formatted tables in results/empirical/.\n")
