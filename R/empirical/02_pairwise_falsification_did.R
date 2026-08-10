# 02_pairwise_falsification_did.R
# Pairwise falsification diagnostic and DID analysis
# "Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach"
#
# Run this script from the repository root after:
#   R/empirical/01_conventional_dif_screening.R

if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop(
    'Install the sandwich package first: install.packages("sandwich")'
  )
}

DATA_FILE <- file.path(
  "data", "derived", "piaac_korea_locator_analysis.csv"
)

SCREENING_FILE <- file.path(
  "results", "empirical", "conventional_dif_screening.csv"
)

if (!file.exists(DATA_FILE)) {
  stop(
    "Prepared data not found. Run ",
    "R/empirical/00_prepare_empirical_data.R first."
  )
}

if (!file.exists(SCREENING_FILE)) {
  stop(
    "Screening results not found. Run ",
    "R/empirical/01_conventional_dif_screening.R first."
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

ITEM_NAMES <- c(
  "C601C06S",
  "C815P001S",
  "C815P002S",
  "C832P001S",
  "C832P002S",
  "C813P001S",
  "C833P001S",
  "C833P002S"
)

ALPHA_LEVEL <- 0.05
CONFIDENCE_LEVEL <- 0.95

candidate_anchors <- screening$Item[
  screening$Candidate_anchor
]

test_items <- setdiff(
  ITEM_NAMES,
  candidate_anchors
)

if (length(candidate_anchors) == 0L) {
  stop("No candidate anchors were retained by the screening stage.")
}

if (length(test_items) == 0L) {
  stop("No test items remain after candidate-anchor screening.")
}


# ------------------------------------------------------------------------------
# HC3 regression for one anchor--test pair
# ------------------------------------------------------------------------------

analyze_pair <- function(
    test_item,
    anchor_item,
    data,
    confidence_level = CONFIDENCE_LEVEL
) {
  pair_data <- data.frame(
    G = as.integer(data$G),
    response_difference =
      data[[test_item]] - data[[anchor_item]]
  )

  fit <- stats::lm(
    response_difference ~ G,
    data = pair_data
  )

  robust_vcov <- sandwich::vcovHC(
    fit,
    type = "HC3"
  )

  estimates <- stats::coef(fit)
  robust_se <- sqrt(diag(robust_vcov))

  residual_df <- stats::df.residual(fit)

  critical_value <- stats::qt(
    1 - (1 - confidence_level) / 2,
    df = residual_df
  )

  robust_p_value <- function(estimate, se) {
    statistic <- estimate / se

    2 * stats::pt(
      abs(statistic),
      df = residual_df,
      lower.tail = FALSE
    )
  }

  alpha_estimate <- unname(estimates["(Intercept)"])
  alpha_se <- unname(robust_se["(Intercept)"])
  alpha_p <- robust_p_value(alpha_estimate, alpha_se)

  did_estimate <- unname(estimates["G"])
  did_se <- unname(robust_se["G"])
  did_p <- robust_p_value(did_estimate, did_se)

  data.frame(
    Test_item = test_item,
    Anchor_item = anchor_item,
    N = nrow(pair_data),
    alpha = alpha_estimate,
    SE_alpha = alpha_se,
    CI_alpha_lower = alpha_estimate - critical_value * alpha_se,
    CI_alpha_upper = alpha_estimate + critical_value * alpha_se,
    p_alpha = alpha_p,
    Diagnostic = ifelse(
      !is.na(alpha_p) & alpha_p >= ALPHA_LEVEL,
      "Did not reject",
      "Rejected"
    ),
    DID = did_estimate,
    SE_DID = did_se,
    CI_lower = did_estimate - critical_value * did_se,
    CI_upper = did_estimate + critical_value * did_se,
    p_DID = did_p,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# Analyze every remaining test item x candidate anchor pair
# ------------------------------------------------------------------------------

pair_results <- do.call(
  rbind,
  lapply(test_items, function(test_item) {
    do.call(
      rbind,
      lapply(candidate_anchors, function(anchor_item) {
        analyze_pair(
          test_item = test_item,
          anchor_item = anchor_item,
          data = dat
        )
      })
    )
  })
)

pair_results$Test_item <- factor(
  pair_results$Test_item,
  levels = ITEM_NAMES
)

pair_results$Anchor_item <- factor(
  pair_results$Anchor_item,
  levels = ITEM_NAMES
)

pair_results <- pair_results[
  order(pair_results$Test_item, pair_results$Anchor_item),
]

pair_results$Test_item <- as.character(pair_results$Test_item)
pair_results$Anchor_item <- as.character(pair_results$Anchor_item)

row.names(pair_results) <- NULL

dir.create(
  file.path("results", "empirical"),
  recursive = TRUE,
  showWarnings = FALSE
)

utils::write.csv(
  pair_results,
  file.path(
    "results", "empirical", "pairwise_falsification_did.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(utils::sessionInfo()),
  con = file.path(
    "results", "empirical", "sessionInfo_pairwise_did.txt"
  )
)

cat("\nPairwise falsification diagnostic and DID results:\n")
print(pair_results, row.names = FALSE)
