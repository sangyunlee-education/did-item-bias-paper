# 01_conventional_dif_screening.R
# Conventional DIF screening for the empirical illustration
# "Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach"
#
# Run this script from the repository root after:
#   R/empirical/00_prepare_empirical_data.R

required_packages <- c(
  "dplyr", "purrr", "tibble", "mirt", "difR"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(mirt)

DATA_FILE <- file.path(
  "data", "derived", "piaac_korea_locator_analysis.csv"
)

if (!file.exists(DATA_FILE)) {
  stop(
    "Prepared data not found. Run ",
    "R/empirical/00_prepare_empirical_data.R first."
  )
}

dat <- utils::read.csv(
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

ALPHA_LEVEL <- 0.05

if (!all(dat$G %in% c(0L, 1L))) {
  stop("G must be coded 0 = reference and 1 = focal.")
}


# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

flag_from_p <- function(p, alpha_level = ALPHA_LEVEL) {
  ifelse(is.na(p), NA, p < alpha_level)
}

extract_mirt_anova_p <- function(anova_table, row = 2L) {
  if (is.null(anova_table) || nrow(anova_table) < row) {
    return(NA_real_)
  }

  preferred_names <- c("p", "Pr(>Chisq)", "Pr(>Chi)", "Pr(>X2)")
  available_name <- preferred_names[preferred_names %in% names(anova_table)]

  if (length(available_name) > 0L) {
    return(as.numeric(anova_table[[available_name[1L]]][row]))
  }

  p_columns <- grep("^p$|^Pr", names(anova_table), ignore.case = TRUE)

  if (length(p_columns) == 0L) {
    return(NA_real_)
  }

  as.numeric(anova_table[[p_columns[1L]]][row])
}


# ------------------------------------------------------------------------------
# 1. Mantel-Haenszel DIF
# ------------------------------------------------------------------------------

run_mh_screen <- function(data, items, group_var = "G") {
  item_data <- data |>
    dplyr::select(dplyr::all_of(items))

  group <- data[[group_var]]

  purrr::map_dfr(items, function(studied_item) {
    matching_items <- setdiff(items, studied_item)
    matching_score <- rowSums(item_data[, matching_items, drop = FALSE])

    contingency_table <- table(
      factor(item_data[[studied_item]], levels = c(0, 1)),
      factor(group, levels = c(0, 1)),
      factor(matching_score, levels = sort(unique(matching_score)))
    )

    # IMPORTANT: correct = FALSE reproduces the manuscript analysis.
    fit <- tryCatch(
      stats::mantelhaen.test(
        contingency_table,
        correct = FALSE
      ),
      error = function(e) {
        warning(
          "Mantel-Haenszel analysis failed for ",
          studied_item, ": ", conditionMessage(e)
        )
        NULL
      }
    )

    p_value <- if (is.null(fit)) {
      NA_real_
    } else {
      as.numeric(fit$p.value)
    }

    tibble::tibble(
      Item = studied_item,
      Mantel_Haenszel = p_value
    )
  })
}


# ------------------------------------------------------------------------------
# 2. Logistic-regression DIF
# ------------------------------------------------------------------------------

run_logistic_screen <- function(data, items, group_var = "G") {
  item_data <- data |>
    dplyr::select(dplyr::all_of(items))

  group <- data[[group_var]]

  purrr::map_dfr(items, function(studied_item) {
    matching_items <- setdiff(items, studied_item)
    matching_score <- rowSums(item_data[, matching_items, drop = FALSE])

    model_data <- tibble::tibble(
      Y = item_data[[studied_item]],
      G = group,
      matching_score_z = as.numeric(scale(matching_score))
    )

    reduced_model <- tryCatch(
      stats::glm(
        Y ~ matching_score_z,
        family = stats::binomial(),
        data = model_data
      ),
      error = function(e) NULL
    )

    full_model <- tryCatch(
      stats::glm(
        Y ~ matching_score_z * G,
        family = stats::binomial(),
        data = model_data
      ),
      error = function(e) NULL
    )

    if (is.null(reduced_model) || is.null(full_model)) {
      warning("Logistic DIF analysis failed for ", studied_item)

      return(
        tibble::tibble(
          Item = studied_item,
          Logistic_DIF = NA_real_
        )
      )
    }

    comparison <- tryCatch(
      stats::anova(
        reduced_model,
        full_model,
        test = "LRT"
      ),
      error = function(e) NULL
    )

    p_value <- if (is.null(comparison)) {
      NA_real_
    } else {
      as.numeric(comparison$`Pr(>Chi)`[2L])
    }

    tibble::tibble(
      Item = studied_item,
      Logistic_DIF = p_value
    )
  })
}


# ------------------------------------------------------------------------------
# 3. SIBTEST
# ------------------------------------------------------------------------------

run_sibtest_screen <- function(data, items, group_var = "G") {
  item_matrix <- as.matrix(
    data |>
      dplyr::select(dplyr::all_of(items))
  )

  group <- data[[group_var]]

  # This reproduces the original analysis: each studied item is evaluated
  # using the remaining items as the matching set.
  fit <- tryCatch(
    difR::sibTest(
      data = item_matrix,
      member = group,
      anchor = seq_along(items),
      type = "udif"
    ),
    error = function(e) {
      warning("SIBTEST failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fit)) {
    return(
      tibble::tibble(
        Item = items,
        SIBTEST = NA_real_
      )
    )
  }

  p_values <- as.numeric(fit$p.value)

  if (length(p_values) != length(items)) {
    stop(
      "SIBTEST returned ", length(p_values),
      " p-values for ", length(items), " items."
    )
  }

  tibble::tibble(
    Item = items,
    SIBTEST = p_values
  )
}


# ------------------------------------------------------------------------------
# 4. 2PL IRT likelihood-ratio DIF
# ------------------------------------------------------------------------------

run_irt_lrt_screen <- function(data, items, group_var = "G") {
  item_data <- data |>
    dplyr::select(dplyr::all_of(items))

  group <- factor(
    data[[group_var]],
    levels = c(0, 1),
    labels = c("reference", "focal")
  )

  # Fully invariant restricted model with group mean and variance free.
  restricted_model <- tryCatch(
    mirt::multipleGroup(
      data = item_data,
      model = 1,
      group = group,
      itemtype = "2PL",
      invariance = c(
        items,
        "free_means",
        "free_var"
      ),
      verbose = FALSE
    ),
    error = function(e) {
      stop(
        "The fully invariant 2PL model failed: ",
        conditionMessage(e)
      )
    }
  )

  purrr::map_dfr(items, function(studied_item) {
    linking_items <- setdiff(items, studied_item)

    # The studied item's 2PL parameters are free across groups;
    # the other seven items remain invariant linking items.
    less_restricted_model <- tryCatch(
      mirt::multipleGroup(
        data = item_data,
        model = 1,
        group = group,
        itemtype = "2PL",
        invariance = c(
          linking_items,
          "free_means",
          "free_var"
        ),
        verbose = FALSE
      ),
      error = function(e) {
        warning(
          "The less restricted 2PL model failed for ",
          studied_item, ": ", conditionMessage(e)
        )
        NULL
      }
    )

    comparison <- if (is.null(less_restricted_model)) {
      NULL
    } else {
      tryCatch(
        anova(
          restricted_model,
          less_restricted_model
        ),
        error = function(e) {
          warning(
            "The 2PL likelihood-ratio comparison failed for ",
            studied_item, ": ", conditionMessage(e)
          )
          NULL
        }
      )
    }

    p_value <- extract_mirt_anova_p(comparison)

    tibble::tibble(
      Item = studied_item,
      IRT_2PL_LRT = p_value
    )
  })
}


# ------------------------------------------------------------------------------
# Run the four DIF screens
# ------------------------------------------------------------------------------

cat("\nRunning Mantel-Haenszel DIF screening...\n")
mh_results <- run_mh_screen(dat, ITEM_NAMES)

cat("Running logistic-regression DIF screening...\n")
logistic_results <- run_logistic_screen(dat, ITEM_NAMES)

cat("Running SIBTEST screening...\n")
sibtest_results <- run_sibtest_screen(dat, ITEM_NAMES)

cat("Running 2PL IRT likelihood-ratio DIF screening...\n")
irt_results <- run_irt_lrt_screen(dat, ITEM_NAMES)

screening <- mh_results |>
  dplyr::full_join(logistic_results, by = "Item") |>
  dplyr::full_join(sibtest_results, by = "Item") |>
  dplyr::full_join(irt_results, by = "Item") |>
  dplyr::mutate(
    all_p_values_available =
      !is.na(Mantel_Haenszel) &
      !is.na(Logistic_DIF) &
      !is.na(SIBTEST) &
      !is.na(IRT_2PL_LRT),
    Candidate_anchor =
      all_p_values_available &
      Mantel_Haenszel >= ALPHA_LEVEL &
      Logistic_DIF >= ALPHA_LEVEL &
      SIBTEST >= ALPHA_LEVEL &
      IRT_2PL_LRT >= ALPHA_LEVEL
  ) |>
  dplyr::arrange(match(Item, ITEM_NAMES))

dir.create(
  file.path("results", "empirical"),
  recursive = TRUE,
  showWarnings = FALSE
)

utils::write.csv(
  screening,
  file.path(
    "results", "empirical", "conventional_dif_screening.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(utils::sessionInfo()),
  con = file.path(
    "results", "empirical", "sessionInfo_screening.txt"
  )
)

cat("\nConventional DIF screening results:\n")
print(screening, row.names = FALSE)

cat("\nCandidate anchors:\n")
print(screening$Item[screening$Candidate_anchor])
