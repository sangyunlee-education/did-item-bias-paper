# 02_pairwise_falsification_did.R
# Pairwise falsification diagnostic and DID analysis
# "Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach"
#
# Run this script from the repository root after
# R/empirical/01_conventional_dif_screening.R.

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

dat <- read.csv(
  DATA_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

screening <- read.csv(
  SCREENING_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

candidate_anchors <- screening$Item[
  screening$Candidate_anchor
]

test_items <- screening$Item[
  !screening$Candidate_anchor
]

if (length(candidate_anchors) == 0) {
  stop("No candidate anchors were retained by the screening stage.")
}

if (length(test_items) == 0) {
  stop("No test items remain after candidate-anchor screening.")
}

ALPHA_LEVEL <- 0.05


# -------------------------------------------------------------------------
# HC3 regression for one anchor--test pair
# -------------------------------------------------------------------------

analyze_pair <- function(test_item, anchor_item, data) {
  pair_data <- data.frame(
    G = as.integer(data$G),
    D = data[[test_item]] - data[[anchor_item]]
  )

  fit <- lm(D ~ G, data = pair_data)

  V_hc3 <- sandwich::vcovHC(
    fit,
    type = "HC3"
  )

  estimates <- coef(fit)
  ses <- sqrt(diag(V_hc3))

  df <- df.residual(fit)
  critical <- qt(0.975, df = df)

  t_stats <- estimates / ses
  p_values <- 2 * pt(
    abs(t_stats),
    df = df,
    lower.tail = FALSE
  )

  ci_lower <- estimates - critical * ses
  ci_upper <- estimates + critical * ses

  # Intercept = alpha = E(Y_T - Y_A | G = 0)
  alpha_hat <- unname(estimates["(Intercept)"])
  se_alpha <- unname(ses["(Intercept)"])
  p_alpha <- unname(p_values["(Intercept)"])

  # Group coefficient = DID
  did_hat <- unname(estimates["G"])
  se_did <- unname(ses["G"])
  p_did <- unname(p_values["G"])

  data.frame(
    Test_item = test_item,
    Anchor_item = anchor_item,
    alpha = alpha_hat,
    SE_alpha = se_alpha,
    p_alpha = p_alpha,
    Diagnostic = ifelse(
      p_alpha < ALPHA_LEVEL,
      "Rejected",
      "Did not reject"
    ),
    DID = did_hat,
    SE_DID = se_did,
    CI_lower = unname(ci_lower["G"]),
    CI_upper = unname(ci_upper["G"]),
    p_DID = p_did,
    stringsAsFactors = FALSE
  )
}


# -------------------------------------------------------------------------
# Analyze every test-item x candidate-anchor pair
# -------------------------------------------------------------------------

pair_grid <- expand.grid(
  Test_item = test_items,
  Anchor_item = candidate_anchors,
  stringsAsFactors = FALSE
)

pair_results <- do.call(
  rbind,
  lapply(seq_len(nrow(pair_grid)), function(i) {
    analyze_pair(
      test_item = pair_grid$Test_item[i],
      anchor_item = pair_grid$Anchor_item[i],
      data = dat
    )
  })
)

# Order results by the test-item order from the screening table, then anchor.
pair_results$Test_item <- factor(
  pair_results$Test_item,
  levels = test_items
)

pair_results$Anchor_item <- factor(
  pair_results$Anchor_item,
  levels = candidate_anchors
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

write.csv(
  pair_results,
  file.path(
    "results", "empirical", "pairwise_falsification_did.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    "results", "empirical", "sessionInfo_pairwise_did.txt"
  )
)

cat("\nPairwise falsification diagnostic and DID results:\n")
print(pair_results, row.names = FALSE)
