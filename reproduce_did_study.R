# Reproduce the two simulation studies and PIAAC illustration (base R >= 3.6).
# Default: update figures and tables from saved CSVs; no new simulations.
# Rscript reproduce_did_study.R presentation
# Rscript reproduce_did_study.R simulations
# Rscript reproduce_did_study.R empirical

# 1. Settings -----------------------------------------------------------------

RUN_MODE <- "presentation" # Update saved results. Use "simulations" for a new run.
N_REP <- 5000L
OUTPUT_DIR <- "results"
PIAAC_FILE <- "prgkorp2.csv"
MAKE_FIGURES <- TRUE
FIGURE_FONT <- "Arial"     # Install this font to preserve the manuscript design.
FINAL_FIGURE_WIDTH_MM <- 144  # Use the same inclusion width for Figures 4 and 5.

# Optional: reuse raw RDS results instead of generating new samples.
STUDY1_SAVED <- ""
STUDY2_SAVED <- ""

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1L) RUN_MODE <- args[1L]
if (length(args) > 1L) stop("Supply at most one run mode.")
stopifnot(RUN_MODE %in% c("presentation", "simulations", "study1", "study2", "outputs",
                        "empirical", "all", "smoke"))
RUN_STUDY1 <- RUN_MODE %in% c("presentation", "simulations", "study1", "outputs", "all", "smoke")
RUN_STUDY2 <- RUN_MODE %in% c("presentation", "simulations", "study2", "outputs", "all", "smoke")
RUN_EMPIRICAL <- RUN_MODE %in% c("empirical", "all")
if (RUN_MODE == "smoke") {
  N_REP <- 20L
  OUTPUT_DIR <- paste0(OUTPUT_DIR, "_smoke")
  STUDY1_SAVED <- STUDY2_SAVED <- ""
}
if (RUN_MODE == "outputs") {
  if (!nzchar(STUDY1_SAVED)) STUDY1_SAVED <- file.path(OUTPUT_DIR, "study1", "simulation.rds")
  if (!nzchar(STUDY2_SAVED)) STUDY2_SAVED <- file.path(OUTPUT_DIR, "study2", "simulation.rds")
}
ALPHA <- .05
MU <- -.5
N_VALUES <- c(500L, 1000L, 2000L)
TAU_VALUES <- c(0, -.05, -.10)
ETA_TRUE_VALUES <- c(0, .05, .10, .15)
ETA_MULTIPLIERS <- c(0, .10, .25, .50, .75, 1, 1.25)
EMPIRICAL_ETA <- .033
Z <- qnorm(1 - ALPHA / 2)
stopifnot(getRversion() >= "3.6.0", N_REP >= 2, N_REP == as.integer(N_REP))
if (RUN_EMPIRICAL && !file.exists(PIAAC_FILE)) stop("PIAAC CSV not found: ", PIAAC_FILE)
if (RUN_MODE != "presentation" && RUN_STUDY1 &&
    nzchar(STUDY1_SAVED) && !file.exists(STUDY1_SAVED)) {
  stop("Saved Study 1 results not found: ", STUDY1_SAVED)
}
if (RUN_MODE != "presentation" && RUN_STUDY2 &&
    nzchar(STUDY2_SAVED) && !file.exists(STUDY2_SAVED)) {
  stop("Saved Study 2 results not found: ", STUDY2_SAVED)
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (RUN_MODE != "presentation") RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# 2. Estimation and saved results --------------------------------------------

# Integrate a function over N(mu, 1). Used for calibration and population checks.
normal_mean <- function(fun, mu = MU) {
  integrate(function(theta) fun(theta) * dnorm(theta, mu, 1),
            -Inf, Inf, subdivisions = 1000L, rel.tol = 1e-10,
            abs.tol = 1e-12)$value
}

# DID = mean(Y_T - Y_A | G = 1) - mean(Y_T - Y_A | G = 0).
# HC3 variance: s1^2/(n1-1) + s0^2/(n0-1), for paired item differences.
fit_did <- function(y_test, y_anchor, group) {
  stopifnot(length(y_test) == length(y_anchor), length(y_test) == length(group),
            all(y_test %in% c(0, 1)), all(y_anchor %in% c(0, 1)),
            all(group %in% c(0, 1)))
  difference <- y_test - y_anchor
  d0 <- difference[group == 0]
  d1 <- difference[group == 1]
  n0 <- length(d0)
  n1 <- length(d1)
  if (min(n0, n1) < 2) stop("HC3 requires at least two people in each group.")
  c(estimate = mean(d1) - mean(d0),
    se = sqrt(var(d1) / (n1 - 1) + var(d0) / (n0 - 1)))
}

# Uniform logistic DIF. Record failed fits; do not redraw their datasets.
fit_logistic <- function(y_test, matching, group) {
  warnings <- character(0)
  fit <- tryCatch(withCallingHandlers(
    glm(y_test ~ matching + group, family = binomial(),
        control = glm.control(maxit = 50L)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  reason <- ""
  estimate <- se <- NA_real_
  if (inherits(fit, "error")) {
    reason <- conditionMessage(fit)
  } else if (!isTRUE(fit$converged) || isTRUE(fit$boundary)) {
    reason <- "Nonconvergence or boundary fit"
  } else {
    coefficients <- tryCatch(coef(summary(fit)), error = function(e) NULL)
    if ("group" %in% rownames(coefficients)) {
      estimate <- unname(coefficients["group", "Estimate"])
      se <- unname(coefficients["group", "Std. Error"])
    }
    if (!is.finite(estimate) || !is.finite(se) || se <= 0) {
      reason <- "Invalid group coefficient or standard error"
    }
  }
  if (nzchar(reason)) estimate <- se <- NA_real_
  list(estimate = estimate, se = se, reason = reason,
       warning = paste(unique(warnings), collapse = " | "))
}

# Evaluate a confidence interval for a specified truth. Study 2 uses truth=tau.
# bound=0 gives the ordinary interval; bound>0 widens both endpoints.
interval_summary <- function(estimate, se, truth, bound = 0) {
  valid <- is.finite(estimate) & is.finite(se) & se >= 0
  n_valid <- sum(valid)
  if (n_valid < 2) stop("Fewer than two valid replications.")
  lower <- estimate[valid] - Z * se[valid] - bound
  upper <- estimate[valid] + Z * se[valid] + bound
  covered <- lower <= truth & truth <= upper
  coverage <- mean(covered)
  interval_length <- upper - lower
  rejection <- mean(lower > 0 | upper < 0)
  data.frame(n_total = length(estimate), n_valid = n_valid,
             n_failed = sum(!valid), coverage = coverage,
             coverage_mcse = sqrt(coverage * (1 - coverage) / n_valid),
             mean_length = mean(interval_length),
             length_mcse = sd(interval_length) / sqrt(n_valid),
             rejection = rejection,
             rejection_mcse = sqrt(rejection * (1 - rejection) / n_valid))
}

# Read either the original repository schema or the revised raw-result schema.
read_saved <- function(path, design, study) {
  saved <- readRDS(path)
  old <- saved$design
  stopifnot(is.data.frame(old), is.data.frame(saved$raw),
            "condition_id" %in% names(old), !anyDuplicated(old$condition_id))
  old <- old[order(old$condition_id), ]
  if (!"rho" %in% names(old) && "rho_M" %in% names(old)) old$rho <- old$rho_M
  fields <- c("condition_id", "N", "tau", "gamma",
              if (study == 1L) c("delta", "rho") else "eta_true")
  stopifnot(all(fields %in% names(old)), nrow(old) == nrow(design))
  for (name in fields) {
    stopifnot(isTRUE(all.equal(old[[name]], design[[name]], tolerance = 1e-8)))
  }
  alpha <- if (!is.null(saved$alpha)) saved$alpha else saved$metadata$alpha
  stopifnot(length(alpha) == 1L, isTRUE(all.equal(alpha, ALPHA)))
  if ("mu" %in% names(old)) stopifnot(all(old$mu == MU))
  if (!is.null(saved$mu)) stopifnot(identical(as.numeric(saved$mu), MU))
  if (study == 2L && "violation" %in% names(old)) {
    stopifnot(all(old$violation == "Smooth"), identical(saved$metadata$anchor_irf,
              "eta + (1 - 2 * eta) * plogis(1.5 * theta)"))
  }
  fields <- c("condition_id", "replication", "did", "did_se",
              if (study == 1L) c("logistic", "logistic_se"))
  stopifnot(all(fields %in% names(saved$raw)))
  raw <- saved$raw[, fields]
  stopifnot(!anyNA(raw$condition_id), !anyNA(raw$replication),
            !anyDuplicated(raw[c("condition_id", "replication")]),
            setequal(raw$condition_id, design$condition_id),
            all(table(raw$condition_id) >= 2L),
            length(unique(as.integer(table(raw$condition_id)))) == 1L)
  if (study == 1L) {
    raw$failure <- if ("failure" %in% names(saved$raw)) saved$raw$failure else
      saved$raw$logistic_failure_reason
    raw$warning <- if ("warning" %in% names(saved$raw)) saved$raw$warning else
      ifelse(saved$raw$logistic_warning, "Warning in original run", "")
    stopifnot(length(raw$failure) == nrow(raw), length(raw$warning) == nrow(raw))
  }
  message("Reusing ", path, " (", nrow(raw) / nrow(design), " replications per condition).")
  raw
}

# Calibrate the constant logit shift to the desired probability-scale effect.
if (RUN_MODE != "presentation") {
  calibration <- data.frame(tau = TAU_VALUES, gamma = 0, tau_achieved = 0)
  for (i in seq_len(nrow(calibration))) {
    target <- calibration$tau[i]
    if (target != 0) {
      calibration$gamma[i] <- uniroot(function(gamma) {
        normal_mean(function(theta) plogis(1.5 * theta + gamma) -
                      plogis(1.5 * theta)) - target
      }, c(-1, 1), tol = 1e-12)$root
    }
    gamma <- calibration$gamma[i]
    calibration$tau_achieved[i] <- normal_mean(function(theta) {
      plogis(1.5 * theta + gamma) - plogis(1.5 * theta)
    })
  }
  stopifnot(max(abs(calibration$tau - calibration$tau_achieved)) < 1e-8)
  write.csv(calibration, file.path(OUTPUT_DIR, "gamma_calibration.csv"), row.names = FALSE)
}

# 3. Study 1 -----------------------------------------------------------------

if (RUN_STUDY1 && RUN_MODE != "presentation") {
  folder <- file.path(OUTPUT_DIR, "study1")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  design1 <- expand.grid(N = N_VALUES, tau = TAU_VALUES,
                        delta = c(0, .25, .50), rho = c(1, .8, .6),
                        KEEP.OUT.ATTRS = FALSE)
  design1$condition_id <- seq_len(nrow(design1))
  design1$gamma <- calibration$gamma[match(design1$tau, calibration$tau)]
  design1$seed <- 2026L + 1009L * design1$condition_id

  if (nzchar(STUDY1_SAVED)) {
    raw1 <- read_saved(STUDY1_SAVED, design1, study = 1L)
  } else {
    runs <- vector("list", nrow(design1))
    for (i in seq_len(nrow(design1))) {
      d <- design1[i, ]
      set.seed(d$seed)
      one <- data.frame(condition_id = d$condition_id, replication = seq_len(N_REP),
                        did = NA_real_, did_se = NA_real_, logistic = NA_real_,
                        logistic_se = NA_real_, failure = "", warning = "",
                        stringsAsFactors = FALSE)
      for (r in seq_len(N_REP)) {
        group <- rbinom(d$N, 1, .5)
        theta <- rnorm(d$N, MU * group, 1)
        y_anchor <- rbinom(d$N, 1, plogis(1.5 * theta))
        y_test <- rbinom(d$N, 1, plogis(1.5 * theta + d$gamma * group))
        error <- if (d$rho == 1) numeric(d$N) else
          rnorm(d$N, 0, sqrt((1 - d$rho) / d$rho))
        matching <- theta + d$delta * group + error
        did <- fit_did(y_test, y_anchor, group)
        logistic <- fit_logistic(y_test, matching, group)
        one$did[r] <- did["estimate"]
        one$did_se[r] <- did["se"]
        one$logistic[r] <- logistic$estimate
        one$logistic_se[r] <- logistic$se
        one$failure[r] <- logistic$reason
        one$warning[r] <- logistic$warning
      }
      runs[[i]] <- one
      message("Study 1: ", i, "/", nrow(design1))
    }
    raw1 <- do.call(rbind, runs)
  }
  stopifnot(all(is.finite(raw1$did)), all(is.finite(raw1$did_se)), all(raw1$did_se >= 0))
  if (!nzchar(STUDY1_SAVED)) {
    saveRDS(list(design = design1, raw = raw1, alpha = ALPHA, mu = MU),
            file.path(folder, "simulation.rds"))
  }
  write.csv(design1, file.path(folder, "design.csv"), row.names = FALSE)

  # Pool raw DID estimates over the nine matching-variable conditions.
  # Pooling is valid here because those conditions do not change DID's inputs.
  data1 <- merge(raw1, design1, by = "condition_id")
  keys <- expand.grid(N = N_VALUES, tau = TAU_VALUES)
  did_rows <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    x <- data1[data1$N == keys$N[i] & data1$tau == keys$tau[i], ]
    did_rows[[i]] <- cbind(keys[i, ],
      data.frame(bias = mean(x$did) - keys$tau[i],
                 rmse = sqrt(mean((x$did - keys$tau[i])^2)),
                 se_sd_ratio = mean(x$did_se) / sd(x$did)),
      interval_summary(x$did, x$did_se, keys$tau[i]))
  }
  did_summary <- do.call(rbind, did_rows)
  logistic_rows <- vector("list", nrow(design1))
  for (i in seq_len(nrow(design1))) {
    x <- raw1[raw1$condition_id == i, ]
    valid <- is.finite(x$logistic) & is.finite(x$logistic_se) & x$logistic_se > 0
    rate <- if (any(valid)) mean(abs(x$logistic[valid] / x$logistic_se[valid]) > Z) else NA_real_
    logistic_rows[[i]] <- cbind(design1[i, ], data.frame(
      n_total = nrow(x), n_valid = sum(valid), n_failed = sum(!valid),
      n_warned = sum(nzchar(x$warning)), rejection = rate,
      rejection_mcse = if (any(valid)) sqrt(rate * (1 - rate) / sum(valid)) else NA_real_))
  }
  logistic_summary <- do.call(rbind, logistic_rows)
  if (any(logistic_summary$n_failed > 0)) {
    warning("Some logistic fits failed; see n_failed and n_valid in logistic_summary.csv.")
  }
  write.csv(did_summary, file.path(folder, "did_summary.csv"), row.names = FALSE)
  write.csv(logistic_summary, file.path(folder, "logistic_summary.csv"), row.names = FALSE)
}

# 4. Study 2 -----------------------------------------------------------------

if (RUN_STUDY2 && RUN_MODE != "presentation") {
  folder <- file.path(OUTPUT_DIR, "study2")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  design2 <- expand.grid(N = N_VALUES, tau = TAU_VALUES, eta_true = ETA_TRUE_VALUES,
                        KEEP.OUT.ATTRS = FALSE)
  design2$condition_id <- seq_len(nrow(design2))
  design2$gamma <- calibration$gamma[match(design2$tau, calibration$tau)]
  design2$seed <- 2027L + 1009L * design2$condition_id

  # Separate the actual DGP discrepancy (eta_true) from the analyst's bound
  # (eta_assumed). Neither is estimated in this simulation.
  normal_factor <- 4 * pnorm(abs(MU) / 2) - 2
  error_per_eta <- 2 * normal_mean(function(theta) plogis(1.5 * theta)) - 1
  design2$id_error <- design2$eta_true * error_per_eta
  design2$population_did <- design2$tau + design2$id_error
  design2$bound_normal <- design2$eta_true * normal_factor
  design2$bound_df <- 2 * design2$eta_true
  stopifnot(all(abs(design2$id_error) <= design2$bound_normal + 1e-10),
            all(design2$bound_normal <= design2$bound_df))

  if (nzchar(STUDY2_SAVED)) {
    raw2 <- read_saved(STUDY2_SAVED, design2, study = 2L)
  } else {
    runs <- vector("list", nrow(design2))
    for (i in seq_len(nrow(design2))) {
      d <- design2[i, ]
      set.seed(d$seed)
      one <- data.frame(condition_id = d$condition_id, replication = seq_len(N_REP),
                        did = NA_real_, did_se = NA_real_)
      for (r in seq_len(N_REP)) {
        group <- rbinom(d$N, 1, .5)
        theta <- rnorm(d$N, MU * group, 1)
        p0 <- plogis(1.5 * theta)
        p_anchor <- d$eta_true + (1 - 2 * d$eta_true) * p0
        y_anchor <- rbinom(d$N, 1, p_anchor)
        y_test <- rbinom(d$N, 1, plogis(1.5 * theta + d$gamma * group))
        did <- fit_did(y_test, y_anchor, group)
        one$did[r] <- did["estimate"]
        one$did_se[r] <- did["se"]
      }
      runs[[i]] <- one
      message("Study 2: ", i, "/", nrow(design2))
    }
    raw2 <- do.call(rbind, runs)
  }
  stopifnot(all(is.finite(raw2$did)), all(is.finite(raw2$did_se)), all(raw2$did_se >= 0))
  if (!nzchar(STUDY2_SAVED)) {
    saveRDS(list(design = design2, raw = raw2, alpha = ALPHA, mu = MU,
                 eta_multipliers = ETA_MULTIPLIERS),
            file.path(folder, "simulation.rds"))
  }
  write.csv(design2, file.path(folder, "design.csv"), row.names = FALSE)

  # Analyze every method and eta choice on the SAME simulated datasets.
  # Reanalysis costs little and needs no additional simulation conditions.
  methods <- c("Unadjusted", "Normal", "Distribution-free")
  rows <- list()
  row_id <- 0L
  diagnostics <- vector("list", nrow(design2))
  for (i in seq_len(nrow(design2))) {
    d <- design2[i, ]
    x <- raw2[raw2$condition_id == d$condition_id, ]
    for (multiplier in ETA_MULTIPLIERS) {
      eta_assumed <- multiplier * d$eta_true
      bounds <- c(0, eta_assumed * normal_factor, 2 * eta_assumed)
      for (j in seq_along(methods)) {
        row_id <- row_id + 1L
        rows[[row_id]] <- cbind(d[, c("condition_id", "N", "tau", "eta_true", "id_error")],
          data.frame(method = methods[j], eta_multiplier = multiplier,
                     eta_assumed = eta_assumed, bound = bounds[j],
                     bound_covers_error = bounds[j] >= abs(d$id_error) - 1e-12,
                     stringsAsFactors = FALSE),
          interval_summary(x$did, x$did_se, truth = d$tau, bound = bounds[j]))
      }
    }
    # Check the ordinary CI for its own estimand, population DID.
    # Its coverage can be nominal even when coverage of tau fails.
    diagnostics[[i]] <- cbind(d,
      data.frame(mean_did = mean(x$did), empirical_error = mean(x$did) - d$tau,
                 empirical_error_mcse = sd(x$did) / sqrt(nrow(x))),
      interval_summary(x$did, x$did_se, truth = d$population_did))
  }
  sensitivity <- do.call(rbind, rows)
  main <- sensitivity[sensitivity$eta_multiplier == 1, ]
  write.csv(main, file.path(folder, "main_summary.csv"), row.names = FALSE)
  write.csv(sensitivity, file.path(folder, "eta_summary.csv"), row.names = FALSE)
  write.csv(do.call(rbind, diagnostics), file.path(folder, "did_diagnostics.csv"),
            row.names = FALSE)

  # Relative total length is descriptive; precision alone does not prove validity.
  width <- reshape(main[, c("condition_id", "N", "tau", "eta_true", "method", "mean_length")],
                   idvar = c("condition_id", "N", "tau", "eta_true"),
                   timevar = "method", direction = "wide")
  width$normal_to_df <- width[["mean_length.Normal"]] / width[["mean_length.Distribution-free"]]
  width$normal_percent_shorter <- 100 * (1 - width$normal_to_df)
  write.csv(width, file.path(folder, "length_comparison.csv"), row.names = FALSE)

  # Population thresholds: below these multipliers the chosen bound is too small
  # for the actual identification error in THIS smooth DGP.
  thresholds <- data.frame(method = c("Normal", "Distribution-free"),
    minimum_multiplier = c(abs(error_per_eta) / normal_factor, abs(error_per_eta) / 2))
  write.csv(thresholds, file.path(folder, "eta_thresholds.csv"), row.names = FALSE)
}

# 5. PIAAC illustration ------------------------------------------------------

if (RUN_EMPIRICAL) {
  folder <- file.path(OUTPUT_DIR, "empirical")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  # These are illustrative, unweighted model-based analyses. They do not account
  # for the complex survey design. No design-based inference is claimed.
  piaac <- read.csv(PIAAC_FILE, sep = ";", dec = ".",
                    na.strings = c(".", ".n", ".v"), check.names = FALSE)
  pv_names <- paste0("PVLIT", 1:10)
  required <- c("AGEG10LFS", "E320004S", "E320003S", pv_names)
  if (!all(required %in% names(piaac))) stop("Required PIAAC columns are missing.")
  keep <- piaac$AGEG10LFS %in% c(2, 4) & !is.na(piaac$E320004S) & !is.na(piaac$E320003S)
  dat <- piaac[keep, required]
  dat$group <- as.integer(dat$AGEG10LFS == 4)
  dat$y_test <- dat$E320004S
  dat$y_anchor <- dat$E320003S
  stopifnot(is.numeric(dat$y_test), is.numeric(dat$y_anchor),
            all(dat$y_test %in% c(0, 1)), all(dat$y_anchor %in% c(0, 1)))
  if (nrow(dat) != 536 || sum(dat$group == 0) != 285 || sum(dat$group == 1) != 251) {
    warning("Analytic sample differs from the manuscript; inspect the input data.")
  }

  # Analyze each plausible value separately, then use MI combining rules.
  pv_results <- data.frame(PV = pv_names, n = 0, estimate = 0, variance = 0,
                           mu = 0, variance_ratio = 0, warning = "", stringsAsFactors = FALSE)
  for (i in seq_along(pv_names)) {
    pv <- dat[[pv_names[i]]]
    stopifnot(is.numeric(pv), !any(is.infinite(pv)))
    use <- !is.na(pv)
    group <- dat$group[use]
    reference <- pv[use][group == 0]
    focal <- pv[use][group == 1]
    stopifnot(length(reference) >= 2, length(focal) >= 2, sd(reference) > 0)
    matching <- (pv[use] - mean(reference)) / sd(reference)
    logistic <- fit_logistic(dat$y_test[use], matching, group)
    if (nzchar(logistic$reason)) stop(pv_names[i], ": ", logistic$reason)
    pv_results$n[i] <- sum(use)
    pv_results$estimate[i] <- logistic$estimate
    pv_results$variance[i] <- logistic$se^2
    pv_results$mu[i] <- (mean(focal) - mean(reference)) / sd(reference)
    pv_results$variance_ratio[i] <- var(focal) / var(reference)
    pv_results$warning[i] <- logistic$warning
  }
  m <- nrow(pv_results)
  logistic_estimate <- mean(pv_results$estimate)
  logistic_se <- sqrt(mean(pv_results$variance) + (1 + 1 / m) * var(pv_results$estimate))
  did <- fit_did(dat$y_test, dat$y_anchor, dat$group)
  estimates <- c(logistic_estimate, unname(did["estimate"]))
  standard_errors <- c(logistic_se, unname(did["se"]))
  empirical <- data.frame(method = c("Logistic-regression DIF", "DID"),
    estimate = estimates, se = standard_errors,
    lower = estimates - Z * standard_errors, upper = estimates + Z * standard_errors,
    p = 2 * pnorm(abs(estimates / standard_errors), lower.tail = FALSE))

  # Average the signed PV differences, then use their absolute value.
  # This plug-in mu is a scenario input, not a value known without uncertainty.
  mu_hat <- mean(pv_results$mu)
  bounds <- c(0, EMPIRICAL_ETA * (4 * pnorm(abs(mu_hat) / 2) - 2), 2 * EMPIRICAL_ETA)
  empirical_sensitivity <- data.frame(
    method = c("Unadjusted", "Normal", "Distribution-free"), eta = EMPIRICAL_ETA,
    mu = c(NA, mu_hat, NA), bound = bounds,
    lower = unname(did["estimate"] - Z * did["se"]) - bounds,
    upper = unname(did["estimate"] + Z * did["se"]) + bounds)
  write.csv(empirical, file.path(folder, "empirical_results.csv"), row.names = FALSE)
  write.csv(empirical_sensitivity, file.path(folder, "sensitivity.csv"), row.names = FALSE)
  write.csv(pv_results, file.path(folder, "pv_results.csv"), row.names = FALSE)
  writeLines(c(paste("Input:", normalizePath(PIAAC_FILE)),
               paste("MD5:", unname(tools::md5sum(PIAAC_FILE))),
               paste("N:", nrow(dat)), paste("Reference:", sum(dat$group == 0)),
               paste("Focal:", sum(dat$group == 1)),
               "Unweighted; survey-design uncertainty not included.",
               "Sensitivity inputs treated as fixed; calibration uncertainty not included."),
             file.path(folder, "analysis_notes.txt"))
  print(empirical, row.names = FALSE)
  print(empirical_sensitivity, row.names = FALSE)
}

# Presentation mode reads the saved summaries only; it never simulates samples.
if (RUN_MODE == "presentation") {
  required <- c("study1/did_summary.csv", "study1/logistic_summary.csv",
                "study2/main_summary.csv", "study2/eta_summary.csv")
  missing <- required[!file.exists(file.path(OUTPUT_DIR, required))]
  if (length(missing)) stop(
    "Presentation mode needs saved CSV summaries under OUTPUT_DIR. Missing: ",
    paste(missing, collapse = ", "),
    ". Copy existing summaries to these paths, use outputs with saved RDS files, ",
    "or explicitly run simulations. No simulations were started.")
  did_summary <- read.csv(file.path(OUTPUT_DIR, required[1]))
  logistic_summary <- read.csv(file.path(OUTPUT_DIR, required[2]))
  main <- read.csv(file.path(OUTPUT_DIR, required[3]))
  sensitivity <- read.csv(file.path(OUTPUT_DIR, required[4]))
  stopifnot(all(c("N", "tau", "bias", "rmse", "se_sd_ratio", "coverage",
                  "rejection", "n_total") %in% names(did_summary)),
            all(c("N", "tau", "delta", "rho", "rejection", "n_total",
                  "n_valid", "n_failed") %in% names(logistic_summary)))
  for (x in list(main, sensitivity)) {
    stopifnot(all(c("N", "tau", "eta_true", "method", "coverage", "mean_length",
                    "rejection", "n_total") %in% names(x)),
              all(is.finite(x$coverage)), all(x$coverage >= 0 & x$coverage <= 1),
              all(is.finite(x$rejection)), all(x$rejection >= 0 & x$rejection <= 1),
              all(is.finite(x$mean_length)), all(x$mean_length >= 0))
  }
  stopifnot("eta_multiplier" %in% names(sensitivity))
  message("Reusing four CSV summaries; no simulations or empirical analyses will run.")
}

# 6. Figures -----------------------------------------------------------------

DEVICE_WIDTH <- 9.4
DEVICE_HEIGHT <- 8.85
FIGURE5_HEIGHT <- 11.8  # Three rows, with unchanged text sizes at 144 mm width.
DEVICE_POINTSIZE <- 12
stopifnot(length(FINAL_FIGURE_WIDTH_MM) == 1L,
          is.finite(FINAL_FIGURE_WIDTH_MM), FINAL_FIGURE_WIDTH_MM > 0)
FIGURE_SCALE <- FINAL_FIGURE_WIDTH_MM / (25.4 * DEVICE_WIDTH)

TEXT_PT <- c(tick = 9, legend = 9, title = 9.5, panel = 10.5, x = 10, y = 10)
TEXT_CEX <- TEXT_PT / (DEVICE_POINTSIZE * FIGURE_SCALE)

LINE_PT <- c(grid = .35, reference = .45, secondary = .8, primary = 1,
             marker = .6, axis = .5)
LINE_LWD <- LINE_PT / (.75 * FIGURE_SCALE)
POINT_CEX <- 1.1

METHODS <- c("Unadjusted", "Normal", "Distribution-free")
COLORS <- c("grey45", "black", "grey30")
LINES <- c(2, 1, 4)
POINTS <- c(21, 21, 24)
FILLS <- c("white", "black", "white")
WIDTHS <- unname(LINE_LWD[c("secondary", "primary", "secondary")])

open_figure <- function(path, height = DEVICE_HEIGHT) {
  if (identical(Sys.info()[["sysname"]], "Darwin") && isTRUE(capabilities("aqua"))) {
    grDevices::quartz(type = "pdf", file = path, width = DEVICE_WIDTH, height = height,
                     family = FIGURE_FONT, pointsize = DEVICE_POINTSIZE)
  } else if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(path, width = DEVICE_WIDTH, height = height,
                        family = FIGURE_FONT, pointsize = DEVICE_POINTSIZE, onefile = TRUE)
  } else {
    stop("PDF output requires native macOS Quartz or Cairo support in R.")
  }
}

blank_strip <- function() {
  par(mar = c(0, 0, 0, 0), cex = 1, xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
}

# Figure 4: rejection rates under invalidity (A) and unreliability (B).
folder <- file.path(OUTPUT_DIR, "study1")
if (MAKE_FIGURES && RUN_STUDY1) {
  did <- read.csv(file.path(folder, "did_summary.csv"))
  logistic <- read.csv(file.path(folder, "logistic_summary.csv"))
  taus <- c(0, -.05, -.10)
  open_figure(file.path(folder, "figure4.pdf"))
  par(oma = c(.4, 4.8, .4, .3), family = FIGURE_FONT,
      ps = DEVICE_POINTSIZE, lwd = LINE_LWD["axis"])
  layout(rbind(c(1, 1, 1), c(2, 3, 4), c(5, 5, 5), c(6, 6, 6),
               c(0, 0, 0),
               c(7, 7, 7), c(8, 9, 10), c(11, 11, 11), c(12, 12, 12)),
         heights = c(.32, 2.3, .24, .48, .10, .32, 2.3, .24, .48))
  colors <- c("black", "grey35", "grey45", "grey55")
  line_types <- c(1, 2, 4, 5)
  symbols <- c(21, 21, 24, 22)
  widths <- unname(LINE_LWD[c("primary", "secondary", "secondary", "secondary")])
  for (block in 1:2) {
    blank_strip()
    text(0, .5, if (block == 1) "A" else "B", adj = c(0, .5), cex = TEXT_CEX["panel"], font = 2)
    for (n in c(500, 1000, 2000)) {
      par(mar = c(2.3, 2.8, 2.5, .8), mgp = c(1.7, .5, 0), tcl = -.2,
          las = 1, xaxs = "i", yaxs = "i", cex = 1)
      plot(NA, xlim = c(.005, -.105), ylim = c(-.025, 1.025),
           xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "l")
      abline(h = c(.25, .5, .75, 1), col = "grey90", lwd = LINE_LWD["grid"])
      abline(h = .05, col = "grey60", lty = 3, lwd = LINE_LWD["reference"])
      axis(1, at = taus, labels = c("0", "-0.05", "-0.10"), cex.axis = TEXT_CEX["tick"])
      axis(2, at = seq(0, 1, .2), labels = c("0", ".2", ".4", ".6", ".8", "1.0"), cex.axis = TEXT_CEX["tick"])
      mtext(bquote(italic(N) == .(format(n, big.mark = ",", trim = TRUE))),
            side = 3, line = .8, cex = TEXT_CEX["title"], las = 1)
      z <- did[did$N == n, ]
      series <- list(z$rejection[match(taus, z$tau)])
      for (k in 1:3) {
        delta <- if (block == 1) c(0, .25, .5)[k] else 0
        rho <- if (block == 1) 1 else c(1, .8, .6)[k]
        z <- logistic[logistic$N == n & logistic$delta == delta & logistic$rho == rho, ]
        series[[k + 1]] <- z$rejection[match(taus, z$tau)]
      }
      for (j in c(2, 3, 4, 1)) {
        lines(taus, series[[j]], col = colors[j], lty = line_types[j], lwd = widths[j])
      }
      for (j in c(2, 3, 4, 1)) {
        points(taus, series[[j]], col = colors[j], pch = symbols[j],
               bg = if (j == 1) "black" else "white", cex = POINT_CEX, lwd = LINE_LWD["marker"])
      }
    }
    blank_strip()
    text(.5, .5, expression(tau), cex = TEXT_CEX["x"])
    blank_strip()
    labels <- if (block == 1) {
      expression(DID, paste("Logistic: ", delta == 0),
                 paste("Logistic: ", delta == .25), paste("Logistic: ", delta == .50))
    } else {
      expression(DID, paste("Logistic: ", rho == 1.00),
                 paste("Logistic: ", rho == .80), paste("Logistic: ", rho == .60))
    }
    legend("center", labels, horiz = TRUE, bty = "n", cex = TEXT_CEX["legend"],
           col = colors, lty = line_types, lwd = widths, pch = symbols,
           pt.bg = c("black", rep("white", 3)), pt.cex = POINT_CEX,
           seg.len = 1.5, x.intersp = .65)
  }
  mtext("Rejection rate", side = 2, outer = TRUE, line = 2.5,
        las = 0, at = .55, cex = TEXT_CEX["y"])
  dev.off()
}

# Figure 5: coverage (A), mean interval length (B), and rejection rates (C).
folder <- file.path(OUTPUT_DIR, "study2")
if (MAKE_FIGURES && RUN_STUDY2) {
  main <- read.csv(file.path(folder, "main_summary.csv"))
  sensitivity <- read.csv(file.path(folder, "eta_summary.csv"))
  # Use common limits across all effect sizes; keep .95 and 1 clearly visible.
  coverage_limits <- c(max(0, floor((min(main$coverage) - .02) * 10) / 10), 1.025)
  length_limits <- c(0, max(main$mean_length) * 1.05)
  open_figure(file.path(folder, "figure5.pdf"), height = FIGURE5_HEIGHT)
  data <- main[main$tau == -.05, ]
  par(oma = c(.4, 4.8, .4, .3), family = FIGURE_FONT,
      ps = DEVICE_POINTSIZE, lwd = LINE_LWD["axis"])
  layout(rbind(c(1, 1, 1), c(2, 3, 4), c(5, 5, 5), c(6, 6, 6),
               c(0, 0, 0),
               c(7, 7, 7), c(8, 9, 10), c(11, 11, 11), c(12, 12, 12),
               c(0, 0, 0),
               c(13, 13, 13), c(14, 15, 16), c(17, 17, 17), c(18, 18, 18)),
         heights = c(.32, 2.3, .24, .48, .10,
                     .32, 2.3, .24, .48, .10,
                     .32, 2.3, .24, .48))
  x_values <- sort(unique(data$eta_true))
  stopifnot(length(x_values) == 4L, all(c("coverage", "mean_length", "rejection") %in% names(data)))
  x_padding <- diff(range(x_values)) / 22
  x_tick_labels <- sprintf("%.2f", x_values)
  x_tick_labels[x_values == 0] <- "0"
  metrics <- c("coverage", "mean_length", "rejection")
  y_labels <- c("Coverage", "Mean interval length", "Rejection rate")
  y_centers <- numeric(3)
  for (block in 1:3) {
    blank_strip()
    text(0, .5, LETTERS[block], adj = c(0, .5),
         cex = TEXT_CEX["panel"], font = 2)
    metric <- metrics[block]
    limits <- list(coverage_limits, length_limits, c(-.025, 1.025))[[block]]
    ticks <- if (block == 1) {
      seq(ceiling(limits[1] * 10) / 10, 1, by = .1)
    } else if (block == 2) {
      pretty(c(0, limits[2]), n = 4)
    } else seq(0, 1, by = .2)
    ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
    labels <- sprintf("%.1f", ticks)
    if (block != 2) labels <- sub("^0\\.", ".", labels)
    labels[ticks == 0] <- "0"
    for (n in N_VALUES) {
      panel <- data[data$N == n, ]
      par(mar = c(2.3, 2.8, 2.5, .8), mgp = c(1.7, .5, 0),
          tcl = -.2, las = 1, xaxs = "i", yaxs = "i", cex = 1)
      plot(NA, xlim = range(x_values) + c(-1, 1) * x_padding,
           ylim = limits, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "l")
      # Outer-margin coordinates, centered on the actual plotting region.
      figure <- par("fig")
      plot_region <- par("plt")
      center <- figure[3] + mean(plot_region[3:4]) * diff(figure[3:4])
      inner <- par("omd")
      y_centers[block] <- (center - inner[3]) / diff(inner[3:4])
      abline(h = ticks[ticks > 0], col = "grey90", lwd = LINE_LWD["grid"])
      if (block == 1) abline(h = .95, col = "grey60", lty = 3, lwd = LINE_LWD["reference"])
      axis(1, at = x_values, labels = x_tick_labels, cex.axis = TEXT_CEX["tick"])
      axis(2, at = ticks, labels = labels, cex.axis = TEXT_CEX["tick"])
      mtext(bquote(italic(N) == .(format(n, big.mark = ",", trim = TRUE))),
            side = 3, line = .8, cex = TEXT_CEX["title"], las = 1)
      for (j in c(1, 3, 2)) {
        z <- panel[panel$method == METHODS[j], ]
        z <- z[order(z$eta_true), ]
        stopifnot(nrow(z) == length(x_values), !anyDuplicated(z$eta_true),
                  isTRUE(all.equal(z$eta_true, x_values)), all(is.finite(z[[metric]])))
        lines(x_values, z[[metric]], col = COLORS[j], lty = LINES[j], lwd = WIDTHS[j])
      }
      for (j in c(1, 3, 2)) {
        z <- panel[panel$method == METHODS[j], ]
        z <- z[order(z$eta_true), ]
        points(x_values, z[[metric]], col = COLORS[j], pch = POINTS[j],
               bg = FILLS[j], cex = POINT_CEX, lwd = LINE_LWD["marker"])
      }
    }
    blank_strip()
    text(.5, .5, expression(eta[0]), cex = TEXT_CEX["x"])
    blank_strip()
    labels <- c("Unadjusted", "Normal-Distribution Bound", "Distribution-Free Bound")
    legend_widths <- vapply(seq_along(labels), function(j) {
      legend(0, .5, labels[j], xjust = 0, yjust = .5, bty = "n",
             cex = TEXT_CEX["legend"], col = COLORS[j], lty = LINES[j],
             lwd = WIDTHS[j], pch = POINTS[j], pt.bg = FILLS[j],
             pt.cex = POINT_CEX, seg.len = 1.5, x.intersp = .65,
             plot = FALSE)$rect$w
    }, numeric(1))
    gap <- .01
    total <- sum(legend_widths) + 2 * gap
    if (total > 1) stop("Legend exceeds available width; check the installed font.")
    left <- (1 - total) / 2
    for (j in seq_along(labels)) {
      legend(left, .5, labels[j], xjust = 0, yjust = .5, bty = "n",
             cex = TEXT_CEX["legend"], col = COLORS[j], lty = LINES[j],
             lwd = WIDTHS[j], pch = POINTS[j], pt.bg = FILLS[j],
             pt.cex = POINT_CEX, seg.len = 1.5, x.intersp = .65)
      left <- left + legend_widths[j] + gap
    }
  }
  for (block in 1:3) {
    mtext(y_labels[block], side = 2, outer = TRUE, line = 2.5,
          las = 0, at = y_centers[block], cex = TEXT_CEX["y"])
  }
  dev.off()

}

# 7. Supplementary tables ----------------------------------------------------

fmt <- function(x, digits = 3L, leading_zero = FALSE) {
  x <- round(x, digits)
  x[x == 0 & !is.na(x)] <- 0
  out <- sprintf(paste0("%.", digits, "f"), x)
  if (!leading_zero) out <- sub("^(-?)0\\.", "\\1.", out)
  out[is.na(x)] <- "---"
  out
}

write_table <- function(data, headers, caption, label, note, path,
                        breaks = integer(0), panel_titles = NULL) {
  body <- character(0)
  for (i in seq_len(nrow(data))) {
    if (i %in% breaks) body <- c(body, "\\addlinespace")
    title <- panel_titles[as.character(i)]
    if (length(title) && !is.na(title)) {
      body <- c(body, paste0("\\multicolumn{", ncol(data), "}{l}{\\textit{", title,
                             "}} \\\\"), "\\addlinespace")
    }
    body <- c(body, paste0(paste(as.character(unlist(data[i, ], use.names = FALSE)),
                                collapse = " & "), " \\\\"))
  }
  writeLines(c("\\begin{table}[!htbp]", "\\tabcolsep=0pt",
    paste0("\\TBL{\\caption{", caption, "\\label{", label, "}}}"),
    "{\\begin{fntable}",
    paste0("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}",
           paste(rep("r", ncol(data)), collapse = ""), "@{}}"),
    "\\toprule", headers, "\\midrule", body, "\\botrule", "\\end{tabular*}",
    paste0("\\footnotetext[]{\\textit{Note:} ", note, "}"),
    "\\end{fntable}}", "\\end{table}"), path)
}

# S1: DID performance pooled over the nine matching-variable conditions.
if (RUN_STUDY1) {
  folder <- file.path(OUTPUT_DIR, "study1")
  d <- did_summary[order(did_summary$N, -did_summary$tau), ]
  tab <- data.frame(N = format(d$N, big.mark = ",", trim = TRUE),
    tau = fmt(d$tau, 2), bias = fmt(d$bias, 4, leading_zero = TRUE), RMSE = fmt(d$rmse, leading_zero = TRUE),
    SE_SD = fmt(d$se_sd_ratio, leading_zero = TRUE), coverage = fmt(d$coverage, 4),
    rejection = fmt(d$rejection, 4))
  tab$N[duplicated(d$N)] <- ""
  pooled <- unique(d$n_total)
  stopifnot(length(pooled) == 1L)
  note <- paste0(
    "$N$ = sample size; $\\tau$ = true average item-bias effect. ",
    "RMSE = root-mean-square error. Bias and RMSE are relative to $\\tau$. ",
    "SE/SD is the ratio of the mean estimated standard error to the empirical ",
    "standard deviation of the pooled estimates. Coverage is the proportion ",
    "of nominal 95\\% confidence intervals containing $\\tau$. ",
    "Rejection is the proportion of replications in which the two-sided test ",
    "rejects at the .05 level, representing the Type~I error rate at ",
    "$\\tau=0$ and power otherwise. Each row pools ",
    format(pooled, big.mark = ",", trim = TRUE),
    " replications across the nine matching-variable conditions, which affect ",
    "neither the DID estimator nor the distribution of its inputs.")
  write_table(tab,
    "$N$ & $\\tau$ & Bias & RMSE & SE/SD & Coverage & Rejection \\\\",
    "DID Performance in Study~1", "tab:sim1_did", note,
    file.path(folder, "table_S1_did.tex"), breaks = c(4L, 7L))

  # S2: the full crossed logistic-regression DIF design.
  keys <- expand.grid(delta = c(0, .25, .50), N = N_VALUES)
  tab <- data.frame(N = format(keys$N, big.mark = ",", trim = TRUE),
                    delta = fmt(keys$delta, 2, leading_zero = TRUE))
  tab$N[duplicated(keys$N)] <- ""
  for (tau in TAU_VALUES) for (rho in c(1, .8, .6)) {
    z <- logistic_summary[logistic_summary$tau == tau & logistic_summary$rho == rho, ]
    index <- match(paste(keys$N, keys$delta), paste(z$N, z$delta))
    stopifnot(!anyNA(index))
    tab[[paste0("rate", ncol(tab))]] <- fmt(z$rejection[index])
  }
  repetitions <- unique(logistic_summary$n_total)
  stopifnot(length(repetitions) == 1L)
  note <- paste0(
    "$N$ = sample size; $\\tau$ = true average item-bias effect; ",
    "$\\delta$ = direct group effect on the matching variable; ",
    "$\\rho$ = matching-variable reliability. ",
    "Each cell reports the proportion of valid replications in which the ",
    "two-sided test rejects at the .05 level, representing the Type~I error ",
    "rate at $\\tau=0$ and power otherwise. Each condition uses ",
    format(repetitions, big.mark = ",", trim = TRUE), " simulated samples.")
  if (any(logistic_summary$n_failed > 0)) {
    note <- paste0(note, " Failed fits are excluded; valid counts range from ",
      min(logistic_summary$n_valid), " to ", max(logistic_summary$n_valid), ".")
  }
  write_table(tab, c(
    " & & \\multicolumn{3}{c}{$\\tau=0$} & \\multicolumn{3}{c}{$\\tau=-.05$} & \\multicolumn{3}{c}{$\\tau=-.10$} \\\\",
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
    "$N$ & $\\delta$ & \\multicolumn{3}{c}{$\\rho$} & \\multicolumn{3}{c}{$\\rho$} & \\multicolumn{3}{c}{$\\rho$} \\\\",
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
    "of intervals containing $\\tau$; length is mean interval length. ",
    "Rejection is the proportion of intervals excluding zero")
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
    tab[[paste0("rejection", ncol(tab))]] <- fmt(z$rejection[index], 4)
  }
  note <- paste0(
    "$N$ = sample size; $\\eta_0$ = supremum of the absolute difference ",
    "between the test- and anchor-item IRFs under the reference-group ",
    "condition; $\\tau$ = true average item-bias effect. ",
    interval_note, ", representing the Type~I error rate at $\\tau=0$ ",
    "and power otherwise. Each condition is based on ", repetitions,
    " replications, with $\\eta=\\eta_0$ and $\\mu=-0.5$. ",
    "The same simulated samples were used for all three intervals.")
  write_table(tab, c(
    " & & \\multicolumn{3}{c}{Unadjusted} & \\multicolumn{3}{c}{Normal-distribution} & \\multicolumn{3}{c}{Distribution-free} \\\\",
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
    "$N$ & $\\eta_0$ & Coverage & Length & Rejection & Coverage & Length & Rejection & Coverage & Length & Rejection \\\\") ,
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
    tab[[paste0("rejection", ncol(tab))]] <- fmt(z$rejection[index], 4)
  }
  note <- paste0(
    "$N$ = sample size; $\\eta/\\eta_0$ = ratio of the specified sensitivity ",
    "parameter to the supremum of the absolute difference between the test- ",
    "and anchor-item IRFs under the reference-group condition. ",
    "$\\tau$ = true average item-bias effect. ",
    interval_note, ", representing power at $\\tau=-.05$. ",
    "For each $N$, the same ", repetitions,
    " simulated samples were used across all three intervals and values of $\\eta$, with ",
    "$\\tau=-.05$, $\\eta_0=.15$, and $\\mu=-0.5$. ",
    "The unadjusted interval does not depend on $\\eta$.")
  write_table(tab, c(
    " & & \\multicolumn{3}{c}{Unadjusted} & \\multicolumn{3}{c}{Normal-distribution} & \\multicolumn{3}{c}{Distribution-free} \\\\",
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
    "$N$ & $\\eta/\\eta_0$ & Coverage & Length & Rejection & Coverage & Length & Rejection & Coverage & Length & Rejection \\\\") ,
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
message("Finished. Results: ", normalizePath(OUTPUT_DIR))
