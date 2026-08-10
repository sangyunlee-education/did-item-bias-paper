# 00_prepare_empirical_data.R
# Data preparation for the empirical illustration
# "Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach"
#
# Run this script from the repository root.
#
# IMPORTANT:
# The manuscript specifies the eight locator-test items and the final coding
# (correct = 1; incorrect/omission = 0), but it does not record the exact
# raw-PUF gender variable name or all raw response codes. Complete the
# CONFIGURATION block below using the OECD PIAAC Cycle 2 Korea PUF/codebook.

# -------------------------------------------------------------------------
# CONFIGURATION: edit these values once for the raw PUF you are using
# -------------------------------------------------------------------------

RAW_DATA_FILE <- file.path("data", "piaac_cycle2_korea.csv")

# Replace with the exact gender/sex variable name in the PUF.
GROUP_VARIABLE <- "REPLACE_WITH_GROUP_VARIABLE"

# Replace with the raw values identifying males and females.
# Character and numeric codes are both supported.
REFERENCE_VALUE <- "REPLACE_WITH_MALE_VALUE"
FOCAL_VALUE <- "REPLACE_WITH_FEMALE_VALUE"

# The manuscript analysis uses these eight Numeracy locator-test items.
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

# Final analysis coding required by the manuscript:
#   correct response             -> 1
#   incorrect response/omission  -> 0
#
# If the raw PUF variables are already coded 0/1, leave these as written.
# Otherwise replace these vectors with the exact raw codes from the codebook.
CORRECT_VALUES <- c(1)
INCORRECT_OR_OMISSION_VALUES <- c(0)

# Any values not included above are set to NA and are excluded by the
# complete-case rule.
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Read raw PUF
# -------------------------------------------------------------------------

if (!file.exists(RAW_DATA_FILE)) {
  stop(
    "Raw data file not found: ", RAW_DATA_FILE,
    "\nPlace the Korea PIAAC Cycle 2 CSV PUF in data/ or edit RAW_DATA_FILE."
  )
}

raw <- read.csv(
  RAW_DATA_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(GROUP_VARIABLE, ITEM_NAMES)
missing_columns <- setdiff(required_columns, names(raw))

if (length(missing_columns) > 0) {
  stop(
    "The following required columns were not found:\n",
    paste(missing_columns, collapse = ", ")
  )
}


# -------------------------------------------------------------------------
# Recode group membership
# -------------------------------------------------------------------------

raw_group <- raw[[GROUP_VARIABLE]]

G <- rep(NA_integer_, length(raw_group))
G[as.character(raw_group) == as.character(REFERENCE_VALUE)] <- 0L
G[as.character(raw_group) == as.character(FOCAL_VALUE)] <- 1L


# -------------------------------------------------------------------------
# Recode item responses
# -------------------------------------------------------------------------

recode_binary_item <- function(x) {
  out <- rep(NA_integer_, length(x))

  out[as.character(x) %in% as.character(CORRECT_VALUES)] <- 1L
  out[
    as.character(x) %in% as.character(INCORRECT_OR_OMISSION_VALUES)
  ] <- 0L

  out
}

items <- as.data.frame(
  lapply(raw[ITEM_NAMES], recode_binary_item),
  check.names = FALSE
)

analysis_data <- data.frame(
  G = G,
  items,
  check.names = FALSE
)


# -------------------------------------------------------------------------
# Complete-case analytic sample
# -------------------------------------------------------------------------

analysis_data <- analysis_data[
  complete.cases(analysis_data[, c("G", ITEM_NAMES)]),
  c("G", ITEM_NAMES)
]

row.names(analysis_data) <- NULL

if (!all(analysis_data$G %in% c(0L, 1L))) {
  stop("G must contain only 0 (reference) and 1 (focal) after recoding.")
}

for (item in ITEM_NAMES) {
  if (!all(analysis_data[[item]] %in% c(0L, 1L))) {
    stop("Item ", item, " is not binary after recoding.")
  }
}


# -------------------------------------------------------------------------
# Save derived analysis data
# -------------------------------------------------------------------------

dir.create(
  file.path("data", "derived"),
  recursive = TRUE,
  showWarnings = FALSE
)

OUTPUT_FILE <- file.path(
  "data", "derived", "piaac_korea_locator_analysis.csv"
)

write.csv(
  analysis_data,
  OUTPUT_FILE,
  row.names = FALSE
)

cat("\nPrepared empirical analysis data\n")
cat("--------------------------------\n")
cat("Total N:    ", nrow(analysis_data), "\n", sep = "")
cat("Reference:  ", sum(analysis_data$G == 0), "\n", sep = "")
cat("Focal:      ", sum(analysis_data$G == 1), "\n", sep = "")
cat("Saved to:   ", OUTPUT_FILE, "\n\n", sep = "")

# Values reported in the manuscript:
EXPECTED_N <- 5780L
EXPECTED_REFERENCE_N <- 2715L
EXPECTED_FOCAL_N <- 3065L

if (
  nrow(analysis_data) != EXPECTED_N ||
  sum(analysis_data$G == 0) != EXPECTED_REFERENCE_N ||
  sum(analysis_data$G == 1) != EXPECTED_FOCAL_N
) {
  warning(
    "The analytic sample counts do not match the manuscript. ",
    "Check the raw group and item-response coding before proceeding."
  )
}
