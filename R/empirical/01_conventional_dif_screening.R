# 01_conventional_dif_screening.R
# Conventional DIF screening for the empirical illustration
# "Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach"
#
# Methods:
#   1. Mantel-Haenszel
#   2. Logistic-regression DIF
#   3. SIBTEST
#   4. 2PL IRT likelihood-ratio test
#
# Each item is treated as the studied item in turn. The remaining seven
# items provide the basis for matching/linking.
#
# Run this script from the repository root after
# R/empirical/00_prepare_empirical_data.R.

required_packages <- c("difR", "mirt")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    "\nFor example: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

DATA_FILE <- file.path(
  "data", "derived", "piaac_korea_locator_analysis.csv"
)

if (!file.exists(DATA_FILE)) {
  stop(
    "Prepared data not found. Run ",
    "R/empirical/00_prepare_empirical_data.R first."
  )
}

dat <- read.csv(
  DATA_FILE,
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

item_data <- dat[, ITEM_NAMES, drop = FALSE]
G <- as.integer(dat$G)

if (!all(G %in% c(0L, 1L))) {
  stop("G must be coded 0 = reference and 1 = focal.")
}

ALPHA_LEVEL <- 0.05


# -------------------------------------------------------------------------
# 2PL IRT restricted model
# -------------------------------------------------------------------------
#
# The restricted model constrains all eight items to be invariant across
# groups while allowing the focal-group latent mean and variance to differ.
# For each studied item, the comparison model frees that item's 2PL
# parameters and keeps the remaining seven items invariant as linking items.

irt_group <- factor(
  ifelse(G == 0L, "reference", "focal"),
  levels = c("reference", "focal")
)

irt_restricted <- mirt::multipleGroup(
  data = item_data,
  model = 1,
  group = irt_group,
  itemtype = "2PL",
  invariance = c("free_means", "free_var", ITEM_NAMES),
  SE = FALSE,
  verbose = FALSE
)


# -------------------------------------------------------------------------
# Screen each item
# -------------------------------------------------------------------------

results <- vector("list", length(ITEM_NAMES))

for (j in seq_along(ITEM_NAMES)) {
  studied_item <- ITEM_NAMES[j]
  other_indices <- setdiff(seq_along(ITEM_NAMES), j)
  other_items <- ITEM_NAMES[other_indices]

  cat(
    sprintf(
      "[%d/%d] Screening %s\n",
      j, length(ITEM_NAMES), studied_item
    )
  )

  # -----------------------------------------------------------------------
  # 1. Mantel-Haenszel
  # -----------------------------------------------------------------------
  # Matching score = sum of the remaining seven items.

  rest_score <- rowSums(item_data[, other_items, drop = FALSE])

  mh_fit <- difR::difMH(
    Data = item_data,
    group = G,
    focal.name = 1,
    match = rest_score,
    MHstat = "MHChisq",
    correct = TRUE,
    exact = FALSE,
    alpha = ALPHA_LEVEL,
    purify = FALSE,
    p.adjust.method = NULL
  )

  p_mh <- as.numeric(mh_fit$p.value[j])

  # -----------------------------------------------------------------------
  # 2. Logistic-regression DIF
  # -----------------------------------------------------------------------
  # Standardizing the rest score does not change the likelihood-ratio test.
  # Reduced model: item ~ matching score
  # Full model:    item ~ matching score * group
  # The LRT jointly tests the group main effect and group-by-score
  # interaction (2 df).

  y <- item_data[[studied_item]]
  rest_score_z <- as.numeric(scale(rest_score))

  logistic_reduced <- glm(
    y ~ rest_score_z,
    family = binomial(link = "logit")
  )

  logistic_full <- glm(
    y ~ rest_score_z * G,
    family = binomial(link = "logit")
  )

  logistic_lrt <- anova(
    logistic_reduced,
    logistic_full,
    test = "LRT"
  )

  p_logistic <- as.numeric(
    logistic_lrt$`Pr(>Chi)`[nrow(logistic_lrt)]
  )

  # -----------------------------------------------------------------------
  # 3. SIBTEST
  # -----------------------------------------------------------------------
  # The remaining seven items are the matching set.

  sib_fit <- difR::sibTest(
    data = as.matrix(item_data),
    member = G,
    anchor = other_indices,
    type = "udif"
  )

  p_sibtest <- as.numeric(sib_fit$p.value[j])

  # -----------------------------------------------------------------------
  # 4. 2PL IRT likelihood-ratio test
  # -----------------------------------------------------------------------
  # The remaining seven items are invariant linking items. The studied
  # item's slope and intercept are free across groups.

  irt_free_studied <- mirt::multipleGroup(
    data = item_data,
    model = 1,
    group = irt_group,
    itemtype = "2PL",
    invariance = c("free_means", "free_var", other_items),
    SE = FALSE,
    verbose = FALSE
  )

  irt_lrt <- anova(
    irt_restricted,
    irt_free_studied
  )

  p_irt_lrt <- as.numeric(
    irt_lrt[nrow(irt_lrt), "p"]
  )

  results[[j]] <- data.frame(
    Item = studied_item,
    Mantel_Haenszel = p_mh,
    Logistic_DIF = p_logistic,
    SIBTEST = p_sibtest,
    IRT_2PL_LRT = p_irt_lrt,
    stringsAsFactors = FALSE
  )
}

screening <- do.call(rbind, results)

screening$Candidate_anchor <- apply(
  screening[
    c(
      "Mantel_Haenszel",
      "Logistic_DIF",
      "SIBTEST",
      "IRT_2PL_LRT"
    )
  ],
  1,
  function(p) all(p >= ALPHA_LEVEL)
)

dir.create(
  file.path("results", "empirical"),
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  screening,
  file.path(
    "results", "empirical", "conventional_dif_screening.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    "results", "empirical", "sessionInfo_screening.txt"
  )
)

cat("\nConventional DIF screening results:\n")
print(screening, row.names = FALSE)

cat("\nCandidate anchors:\n")
print(screening$Item[screening$Candidate_anchor])
