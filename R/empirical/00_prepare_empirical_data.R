# 00_prepare_empirical_data.R
# Data preparation for the empirical illustration
# "Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach"
#
# Run this script from the repository root.

required_packages <- c("dplyr", "tidyr")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

# Preferred GitHub/repository location for the raw Korea PIAAC Cycle 2 file.
# For convenience, the script also accepts prgkorp2.csv in the repository root.
RAW_DATA_CANDIDATES <- c(
  file.path("data", "piaac", "prgkorp2.csv"),
  "prgkorp2.csv"
)

existing_data_files <- RAW_DATA_CANDIDATES[file.exists(RAW_DATA_CANDIDATES)]

if (length(existing_data_files) == 0L) {
  stop(
    "PIAAC data file not found.\n",
    "Place prgkorp2.csv at data/piaac/prgkorp2.csv ",
    "or in the repository root."
  )
}

RAW_DATA_FILE <- existing_data_files[1L]
DELIMITER <- ";"

GROUP_VARIABLE <- "GENDER_R"
REFERENCE_VALUE <- 1
FOCAL_VALUE <- 2

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

# Raw PIAAC score coding used in the analysis.
CORRECT_VALUES <- "1"
INCORRECT_OR_OMISSION_VALUES <- c("0", "7")

# Manuscript sample counts.
EXPECTED_N <- 5780L
EXPECTED_REFERENCE_N <- 2715L
EXPECTED_FOCAL_N <- 3065L


# ------------------------------------------------------------------------------
# Read and recode the raw PIAAC file
# ------------------------------------------------------------------------------

raw_data <- utils::read.csv(
  RAW_DATA_FILE,
  sep = DELIMITER,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_variables <- c(GROUP_VARIABLE, ITEM_NAMES)
missing_variables <- setdiff(required_variables, names(raw_data))

if (length(missing_variables) > 0L) {
  stop(
    "The following required variables are missing: ",
    paste(missing_variables, collapse = ", ")
  )
}

analysis_data <- raw_data |>
  dplyr::select(
    dplyr::all_of(GROUP_VARIABLE),
    dplyr::all_of(ITEM_NAMES)
  ) |>
  dplyr::mutate(
    G = dplyr::case_when(
      .data[[GROUP_VARIABLE]] == REFERENCE_VALUE ~ 0L,
      .data[[GROUP_VARIABLE]] == FOCAL_VALUE ~ 1L,
      TRUE ~ NA_integer_
    )
  ) |>
  dplyr::select(-dplyr::all_of(GROUP_VARIABLE)) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(ITEM_NAMES),
      ~ dplyr::case_when(
        as.character(.x) %in% CORRECT_VALUES ~ 1,
        as.character(.x) %in% INCORRECT_OR_OMISSION_VALUES ~ 0,
        TRUE ~ NA_real_
      )
    )
  ) |>
  tidyr::drop_na()

observed_groups <- sort(unique(analysis_data$G))

if (!identical(observed_groups, c(0L, 1L))) {
  stop("Both reference (G = 0) and focal (G = 1) groups must be present.")
}


# ------------------------------------------------------------------------------
# Check the analytic sample against the manuscript
# ------------------------------------------------------------------------------

n_total <- nrow(analysis_data)
n_reference <- sum(analysis_data$G == 0L)
n_focal <- sum(analysis_data$G == 1L)

cat("\nPrepared empirical analysis data\n")
cat("--------------------------------\n")
cat("Raw file:   ", RAW_DATA_FILE, "\n", sep = "")
cat("Total N:    ", n_total, "\n", sep = "")
cat("Reference:  ", n_reference, "\n", sep = "")
cat("Focal:      ", n_focal, "\n", sep = "")

if (
  n_total != EXPECTED_N ||
  n_reference != EXPECTED_REFERENCE_N ||
  n_focal != EXPECTED_FOCAL_N
) {
  stop(
    "\nThe analytic sample does not match the manuscript.\n",
    "Expected N = 5780, reference = 2715, focal = 3065.\n",
    "Check that the correct Korea PIAAC Cycle 2 file is being used."
  )
}

cat("Sample check: PASS\n\n")


# ------------------------------------------------------------------------------
# Save derived analysis data
# ------------------------------------------------------------------------------

dir.create(
  file.path("data", "derived"),
  recursive = TRUE,
  showWarnings = FALSE
)

OUTPUT_FILE <- file.path(
  "data", "derived", "piaac_korea_locator_analysis.csv"
)

utils::write.csv(
  analysis_data,
  OUTPUT_FILE,
  row.names = FALSE
)

cat("Saved to: ", OUTPUT_FILE, "\n", sep = "")
