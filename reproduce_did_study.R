    " & & 1.00 & .80 & .60 & 1.00 & .80 & .60 & 1.00 & .80 & .60 \\\\") ,
    "Rejection Rates for Logistic-Regression DIF in Study~1", "tab:sim1_logistic",
    note, file.path(folder, "table_S2_logistic.tex"), breaks = c(4L, 7L))
}

if (RUN_STUDY2) {
  folder <- file.path(OUTPUT_DIR, "study2")
  interval_note <- paste0(
    "Unadjusted denotes the nominal 95\\% DID confidence interval; ",
    "normal-distribution and distribution-free denote the intervals based ",
    "on the corresponding sensitivity bounds. Coverage is the proportion ",
    "of intervals containing $\\tau$; length is mean interval length.")
  repetitions <- unique(main$n_total)
  stopifnot(length(repetitions) == 1L)
  repetitions <- format(repetitions, big.mark = ",", trim = TRUE)

  # S3: additional effect sizes, with eta = eta_0.
  keys <- expand.grid(eta_true = ETA_TRUE_VALUES, N = N_VALUES, tau = c(0, -.10))
  tab <- data.frame(N = format(keys$N, big.mark = ",", trim = TRUE),
                    eta = fmt(keys$eta_true, 2))
  groups <- paste(keys$tau, keys$N)
  tab$N[duplicated(groups)] <- ""
  for (method in METHODS) {
    z <- main[main$method == method, ]
    index <- match(paste(keys$tau, keys$N, keys$eta_true),
                   paste(z$tau, z$N, z$eta_true))
    stopifnot(!anyNA(index))
    tab[[paste0("coverage", ncol(tab))]] <- fmt(z$coverage[index], 4)
    tab[[paste0("length", ncol(tab))]] <- fmt(z$mean_length[index], leading_zero = TRUE)
  }
  note <- paste0(
    "$N$ = sample size; $\\eta_0$ = supremum of the absolute difference ",
    "between the test- and anchor-item IRFs under the reference-group ",
    "condition; $\\tau$ = true average item-bias effect. ",
    interval_note, " Each condition is based on ", repetitions,
    " replications, with $\\eta=\\eta_0$ and $\\mu=-0.5$. ",
    "The same simulated samples were used for all three intervals.")
  write_table(tab, c(
    " & & \\multicolumn{2}{c}{Unadjusted} & \\multicolumn{2}{c}{Normal-distribution} & \\multicolumn{2}{c}{Distribution-free} \\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
    "$N$ & $\\eta_0$ & Coverage & Length & Coverage & Length & Coverage & Length \\\\") ,
    "Interval Performance at Additional Effect Sizes in Study~2",
    "tab:sim2_other_effects", note, file.path(folder, "table_S3_other_effects.tex"),
    breaks = seq(5L, 21L, by = 4L),
    panel_titles = c("1" = "Panel A: $\\tau=0$", "13" = "Panel B: $\\tau=-.10$"))

  # S4: eta choices applied to the same draws at tau = -.05 and eta_0 = .15.
  repetitions <- unique(sensitivity$n_total[sensitivity$tau == -.05 &
                                            sensitivity$eta_true == .15])
  stopifnot(length(repetitions) == 1L)
  repetitions <- format(repetitions, big.mark = ",", trim = TRUE)
  keys <- expand.grid(eta_multiplier = ETA_MULTIPLIERS, N = N_VALUES)
  tab <- data.frame(N = format(keys$N, big.mark = ",", trim = TRUE),
                    ratio = fmt(keys$eta_multiplier, 2, leading_zero = TRUE))
  tab$N[duplicated(keys$N)] <- ""
  for (method in METHODS) {
    z <- sensitivity[sensitivity$tau == -.05 & sensitivity$eta_true == .15 &
                       sensitivity$method == method, ]
    index <- match(paste(keys$N, keys$eta_multiplier), paste(z$N, z$eta_multiplier))
    stopifnot(!anyNA(index))
    tab[[paste0("coverage", ncol(tab))]] <- fmt(z$coverage[index], 4)
    tab[[paste0("length", ncol(tab))]] <- fmt(z$mean_length[index], leading_zero = TRUE)
  }
  note <- paste0(
    "$N$ = sample size; $\\eta/\\eta_0$ = ratio of the specified sensitivity ",
    "parameter to the supremum of the absolute difference between the test- ",
    "and anchor-item IRFs under the reference-group condition. ",
    "$\\tau$ = true average item-bias effect. ",
    interval_note, " ",
    "For each $N$, the same ", repetitions,
    " simulated samples were used across all three intervals and values of $\\eta$, with ",
    "$\\tau=-.05$, $\\eta_0=.15$, and $\\mu=-0.5$. ",
    "The unadjusted interval does not depend on $\\eta$.")
  write_table(tab, c(
    " & & \\multicolumn{2}{c}{Unadjusted} & \\multicolumn{2}{c}{Normal-distribution} & \\multicolumn{2}{c}{Distribution-free} \\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
    "$N$ & $\\eta/\\eta_0$ & Coverage & Length & Coverage & Length & Coverage & Length \\\\") ,
    "Sensitivity to the Specified Value of $\\eta$ in Study~2",
    "tab:sim2_eta_choice", note, file.path(folder, "table_S4_eta_choice.tex"),
    breaks = c(8L, 15L))
}

# 8. Checks and run record ---------------------------------------------------
# Check key identities on the actual results, without generating new datasets.
if (RUN_STUDY2 && RUN_MODE != "presentation") {
  for (i in unique(main$condition_id)) {
    x <- main[main$condition_id == i, ]
    x <- x[match(METHODS, x$method), ]
    stopifnot(all(diff(x$coverage) >= 0),
              max(abs(x$mean_length - x$mean_length[1] - 2 * x$bound)) < 1e-10)
  }
  for (i in seq_len(nrow(design2))) {
    d <- design2[i, ]
    anchor <- function(theta) d$eta_true + (1 - 2 * d$eta_true) * plogis(1.5 * theta)
    population <- normal_mean(function(theta) plogis(1.5 * theta + d$gamma) - anchor(theta)) -
      normal_mean(function(theta) plogis(1.5 * theta) - anchor(theta), mu = 0)
    stopifnot(abs(population - d$population_did) < 1e-8)
  }
}

if (RUN_MODE == "presentation") {
  # Keep the simulation run record intact. These settings describe rendering,
  # not a new analysis or the provenance of the supplied summaries.
  paths <- file.path(OUTPUT_DIR, required)
  writeLines(capture.output(sessionInfo()),
             file.path(OUTPUT_DIR, "presentation_sessionInfo.txt"))
  writeLines(capture.output(dput(list(
    run_mode = RUN_MODE, source_files = required,
    source_md5 = setNames(unname(tools::md5sum(paths)), required),
    figures_written = MAKE_FIGURES, figure_font = FIGURE_FONT,
    final_figure_width_mm = FINAL_FIGURE_WIDTH_MM,
    figure5_height_inches = FIGURE5_HEIGHT,
    run_time = format(Sys.time(), tz = "UTC")))),
    file.path(OUTPUT_DIR, "presentation_settings.txt"))
} else {
  writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
  writeLines(capture.output(dput(list(
  run_mode = RUN_MODE, alpha = ALPHA, mu = MU, n_rep_requested = N_REP,
  figure_font = FIGURE_FONT, final_figure_width_mm = FINAL_FIGURE_WIDTH_MM,
  figure5_height_inches = FIGURE5_HEIGHT,
  N = N_VALUES, tau = TAU_VALUES, eta_true = ETA_TRUE_VALUES,
  eta_multipliers = ETA_MULTIPLIERS, study1_saved = STUDY1_SAVED,
  study2_saved = STUDY2_SAVED, empirical_eta = EMPIRICAL_ETA,
  group_probability = .5, item_responses_conditionally_independent = TRUE,
  intervals_clipped = FALSE, sensitivity_inputs_fixed = TRUE,
  study1_seed = 2026, study2_seed = 2027, seed_increment = 1009,
  run_time = format(Sys.time(), tz = "UTC")))), file.path(OUTPUT_DIR, "settings.txt"))
}
