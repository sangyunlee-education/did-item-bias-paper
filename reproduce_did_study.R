# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach
# Reproduces Studies 1-2, the PIAAC illustration, Figures 4-5, and Tables S1-S4.
# Requires R >= 3.6; no additional R packages are needed.
#
# Set the working directory to this file's folder and provide the Korean PIAAC
# public-use CSV at PIAAC_FILE below. Run the entire file using Source in RStudio
# or Rscript reproduce_did_study.R. Sections follow the manuscript's order.
# All simulations use 5,000 replications per condition and fixed random seeds.
# Outputs are saved under results/study1, results/study2, and results/empirical.
# Figures use Arial and a Quartz or Cairo PDF device; LaTeX tables use CUP macros.

# 1. Inputs and manuscript settings -----------------------------------------
PIAAC_FILE <- "prgkorp2.csv"
OUTPUT_DIR <- "results"
N_REP <- 5000L

# Display settings; detailed figure styling appears in Section 6.
FIGURE_FONT <- "Arial"
FINAL_FIGURE_WIDTH_MM <- 144

# Simulation design: N is sample size, tau is the average item-bias effect,
# and MU is the actual focal-group mean on the reference-group standardized scale.
ALPHA <- .05
MU <- -.5                    # Focal-group latent-trait mean; the reference-group mean is 0.
N_VALUES <- c(500L, 1000L, 2000L)
TAU_VALUES <- c(0, -.05, -.10)
ETA_TRUE_VALUES <- c(.05, .10, .15)
# Supplementary sensitivity inputs, applied to the same samples at tau = -.05
# and actual IRF discrepancy .15. Main analyses use the data-generating values.
ETA_SPECIFIED <- c(.075, .150, .225)
MU_SPECIFIED <- c(-.25, -.50, -.75)
SENSITIVITY_TAU <- -.05
SENSITIVITY_TRUE_ETA <- .15
# Empirical inputs chosen slightly above the reference values .0329 and .72.
EMPIRICAL_ETA <- .035
EMPIRICAL_MU <- -.75
# Internal interval keys; figures and tables set their own display order.
METHODS <- c("Unadjusted", "Normal", "Distribution-free")
Z <- qnorm(1 - ALPHA / 2)

# Check the empirical input before starting the simulations.
stopifnot(getRversion() >= "3.6.0", N_REP >= 2, N_REP == as.integer(N_REP))
if (!file.exists(PIAAC_FILE)) stop("Korean PIAAC CSV not found: ", PIAAC_FILE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 2. Shared statistical calculations ----------------------------------------
# Integrate a function over N(mu, 1). Used for calibration and population checks.
normal_mean <- function(fun, mu = MU) {
  integrate(function(theta) fun(theta) * dnorm(theta, mu, 1),
            -Inf, Inf, subdivisions = 1000L, rel.tol = 1e-10,
            abs.tol = 1e-12)$value
}

# Calibrate the constant logit shift to the desired probability-scale effect.
calibrate_gamma <- function() {
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
  calibration
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

# Estimate the group coefficient in uniform logistic-regression DIF.
# Failed fits are recorded and excluded from its rejection-rate denominator.
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

# Summarize intervals across replications: coverage, mean length, and rejection.
# A positive bound expands both endpoints; rejection means the interval excludes 0.
# The summaries also include Monte Carlo standard errors (mcse).
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

# Format rounded numbers for manuscript tables and the empirical report.
fmt <- function(x, digits = 3L) {
  x <- round(x, digits)
  x[x == 0 & !is.na(x)] <- 0
  out <- sprintf(paste0("%.", digits, "f"), x)
  out <- sub("^(-?)0\\.", "\\1.", out)
  out[is.na(x)] <- "---"
  out
}

# Figures and tables read the CSV summaries produced by the analyses below.
read_summary <- function(relative_path) {
  read.csv(file.path(OUTPUT_DIR, relative_path), stringsAsFactors = FALSE)
}

# 3. Study 1: manuscript Section 4.1 -----------------------------------------

# Vary sample size, item bias, matching-variable invalidity, and reliability.
# The anchor item is valid and IRF equivalence holds.
run_study1 <- function() {
  # Define the crossed design and calibrate the test-item logit shift.
  calibration <- calibrate_gamma()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  folder <- file.path(OUTPUT_DIR, "study1")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  design1 <- expand.grid(N = N_VALUES, tau = TAU_VALUES,
                        delta = c(0, .25, .50), rho = c(1, .8, .6),
                        KEEP.OUT.ATTRS = FALSE)
  design1$condition_id <- seq_len(nrow(design1))
  design1$gamma <- calibration$gamma[match(design1$tau, calibration$tau)]
  design1$seed <- 2026L + 1009L * design1$condition_id

  # Obtain replication-level DID and logistic-regression estimates.
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
  stopifnot(all(is.finite(raw1$did)), all(is.finite(raw1$did_se)), all(raw1$did_se >= 0))
  saveRDS(list(design = design1, raw = raw1, alpha = ALPHA, mu = MU),
          file.path(folder, "simulation.rds"))
  write.csv(calibration, file.path(folder, "gamma_calibration.csv"), row.names = FALSE)
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
  # Summarize logistic-regression rejection rates for each crossed condition.
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

run_study1()

# 4. Study 2: manuscript Section 4.2 -----------------------------------------

# Vary sample size, item bias, and the actual IRF discrepancy.
# Compare unadjusted intervals with two sensitivity-adjusted intervals.
run_study2 <- function() {
  # Define the crossed design and its population identification error.
  calibration <- calibrate_gamma()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  folder <- file.path(OUTPUT_DIR, "study2")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  design2 <- expand.grid(N = N_VALUES, tau = TAU_VALUES, eta_true = ETA_TRUE_VALUES,
                        KEEP.OUT.ATTRS = FALSE)
  # Fixed condition identifiers determine the random seeds.
  design2$condition_id <- 9L + seq_len(nrow(design2))
  design2$gamma <- calibration$gamma[match(design2$tau, calibration$tau)]
  design2$seed <- 2027L + 1009L * design2$condition_id

  # eta_true is the actual maximum IRF difference; eta_assumed is the
  # sensitivity parameter specified for an interval. Both are fixed inputs.
  normal_factor <- 4 * pnorm(abs(MU) / 2) - 2
  error_per_eta <- 2 * normal_mean(function(theta) plogis(1.5 * theta)) - 1
  design2$id_error <- design2$eta_true * error_per_eta
  design2$population_did <- design2$tau + design2$id_error
  design2$bound_normal <- design2$eta_true * normal_factor
  design2$bound_df <- 2 * design2$eta_true
  stopifnot(all(abs(design2$id_error) <= design2$bound_normal + 1e-10),
            all(design2$bound_normal <= design2$bound_df))

  # Obtain DID estimates once; all interval comparisons reuse these samples.
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
  stopifnot(all(is.finite(raw2$did)), all(is.finite(raw2$did_se)), all(raw2$did_se >= 0))
  saveRDS(list(design = design2, raw = raw2, alpha = ALPHA, mu = MU),
          file.path(folder, "simulation.rds"))
  write.csv(calibration, file.path(folder, "gamma_calibration.csv"), row.names = FALSE)
  write.csv(design2, file.path(folder, "design.csv"), row.names = FALSE)

  # Main analysis: the specified eta and mu equal their data-generating values.
  rows <- list()
  row_id <- 0L
  diagnostics <- vector("list", nrow(design2))
  for (i in seq_len(nrow(design2))) {
    d <- design2[i, ]
    x <- raw2[raw2$condition_id == d$condition_id, ]
    bounds <- c(0, d$eta_true * normal_factor, 2 * d$eta_true)
    for (j in seq_along(METHODS)) {
      row_id <- row_id + 1L
      rows[[row_id]] <- cbind(d[, c("condition_id", "N", "tau", "eta_true", "id_error")],
        data.frame(method = METHODS[j], eta_multiplier = 1,
                   eta_assumed = d$eta_true, bound = bounds[j],
                   bound_covers_error = bounds[j] >= abs(d$id_error) - 1e-12,
                   stringsAsFactors = FALSE),
        interval_summary(x$did, x$did_se, truth = d$tau, bound = bounds[j]))
    }
    # Check the ordinary CI for its own estimand, population DID.
    # Its coverage can be nominal even when coverage of tau fails.
    diagnostics[[i]] <- cbind(d,
      data.frame(mean_did = mean(x$did), empirical_error = mean(x$did) - d$tau,
                 empirical_error_mcse = sd(x$did) / sqrt(nrow(x))),
      interval_summary(x$did, x$did_se, truth = d$population_did))
  }
  main <- do.call(rbind, rows)
  write.csv(main, file.path(folder, "main_summary.csv"), row.names = FALSE)
  write.csv(do.call(rbind, diagnostics), file.path(folder, "did_diagnostics.csv"),
            row.names = FALSE)

  # Compare mean lengths of the normal-distribution and distribution-free intervals.
  width <- reshape(main[, c("condition_id", "N", "tau", "eta_true", "method", "mean_length")],
                   idvar = c("condition_id", "N", "tau", "eta_true"),
                   timevar = "method", direction = "wide")
  width$normal_to_df <- width[["mean_length.Normal"]] / width[["mean_length.Distribution-free"]]
  width$normal_percent_shorter <- 100 * (1 - width$normal_to_df)
  write.csv(width, file.path(folder, "length_comparison.csv"), row.names = FALSE)

  # Diagnostic: minimum ratio of specified eta to the actual IRF gap needed
  # for each bound to cover the identification error in this data-generating process.
  thresholds <- data.frame(method = c("Normal", "Distribution-free"),
    minimum_multiplier = c(abs(error_per_eta) / normal_factor, abs(error_per_eta) / 2))
  write.csv(thresholds, file.path(folder, "eta_thresholds.csv"), row.names = FALSE)

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

  # Supplementary analysis: reuse samples at tau = -.05 and actual IRF gap .15.
  # Cross three specified eta values with three specified mu values.
  selected <- design2[abs(design2$tau - SENSITIVITY_TAU) < 1e-10 &
                       abs(design2$eta_true - SENSITIVITY_TRUE_ETA) < 1e-10, ]
  stopifnot(nrow(selected) == length(N_VALUES), setequal(selected$N, N_VALUES),
            !anyDuplicated(selected$N))
  input_grid <- expand.grid(mu_assumed = MU_SPECIFIED, eta_assumed = ETA_SPECIFIED,
                            KEEP.OUT.ATTRS = FALSE)
  rows <- list()
  row_id <- 0L
  for (i in seq_len(nrow(selected))) {
    d <- selected[i, ]
    x <- raw2[raw2$condition_id == d$condition_id, ]
    stopifnot(all(is.finite(x$did)), all(is.finite(x$did_se)), all(x$did_se >= 0))
    for (k in seq_len(nrow(input_grid))) {
      eta <- input_grid$eta_assumed[k]
      mu <- input_grid$mu_assumed[k]
      bounds <- c(0, eta * (4 * pnorm(abs(mu) / 2) - 2), 2 * eta)
      for (j in seq_along(METHODS)) {
        row_id <- row_id + 1L
        rows[[row_id]] <- cbind(
          data.frame(condition_id = d$condition_id, N = d$N, tau = d$tau,
                     eta_true = d$eta_true, mu_true = MU, method = METHODS[j],
                     eta_assumed = eta, mu_assumed = mu, bound = bounds[j],
                     stringsAsFactors = FALSE),
          interval_summary(x$did, x$did_se, truth = d$tau, bound = bounds[j]))
      }
    }
  }
  joint <- do.call(rbind, rows)
  # The CSV includes all settings; S4 displays each distinct interval only once.
  write.csv(joint, file.path(folder, "eta_mu_summary.csv"), row.names = FALSE)
  writeLines(capture.output(dput(list(
    tau = SENSITIVITY_TAU, actual_irf_gap = SENSITIVITY_TRUE_ETA, true_mu = MU,
    eta_specified = ETA_SPECIFIED, mu_specified = MU_SPECIFIED,
    n_rep = unique(joint$n_total),
    source = "new", input_rds = file.path(folder, "simulation.rds")))),
    file.path(folder, "sensitivity_settings.txt"))
  message("Study 2: 27 conditions; nine sensitivity settings on the same samples.")
}

run_study2()

# 5. PIAAC illustration: manuscript Section 5 --------------------------------

run_empirical <- function() {
  if (!file.exists(PIAAC_FILE)) stop("PIAAC CSV not found: ", PIAAC_FILE)
  folder <- file.path(OUTPUT_DIR, "empirical")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  # The illustration is unweighted; standard errors omit the complex survey design.
  piaac <- read.csv(PIAAC_FILE, sep = ";", dec = ".",
                    na.strings = c(".", ".n", ".v"), check.names = FALSE)
  pv_names <- paste0("PVLIT", 1:10)
  required <- c("AGEG10LFS", "E320004S", "E320003S", pv_names)
  if (!all(required %in% names(piaac))) stop("Required PIAAC columns are missing.")
  # Select the two age groups and recode the scored item responses.
  # Raw item scores: 1 = correct, 7 = incorrect, 0 = omitted.
  item_names <- c("E320004S", "E320003S")
  for (item in item_names) {
    scores <- piaac[[item]]
    observed <- scores[!is.na(scores)]
    if (!all(observed %in% c(0, 1, 7))) {
      stop("Unexpected response codes in ", item, "; check the codebook.")
    }
  }

  keep <- piaac$AGEG10LFS %in% c(2, 4) &
    !is.na(piaac$E320004S) &
    !is.na(piaac$E320003S)

  dat <- piaac[keep, required]

  # Reference group: 25-34; focal group: 45-54.
  dat$group <- as.integer(dat$AGEG10LFS == 4)

  # Code correct responses as 1; incorrect and omitted responses as 0.
  dat$y_test <- as.integer(dat$E320004S == 1)
  dat$y_anchor <- as.integer(dat$E320003S == 1)

  if (nrow(dat) != 536 ||
      sum(dat$group == 0) != 285 ||
      sum(dat$group == 1) != 251) {
    stop("Analytic sample differs from the manuscript; inspect the input data.")
  }
  # Fit DIF using each plausible value on its original scale.
  # Pool the ten group coefficients using multiple-imputation combining rules.
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
    matching <- pv[use]
    logistic <- fit_logistic(dat$y_test[use], matching, group)
    if (nzchar(logistic$reason)) stop(pv_names[i], ": ", logistic$reason)
    pv_results$n[i] <- sum(use)
    pv_results$estimate[i] <- logistic$estimate
    pv_results$variance[i] <- logistic$se^2
    pv_results$mu[i] <- (mean(focal) - mean(reference)) / sd(reference)
    pv_results$variance_ratio[i] <- var(focal) / var(reference)
    pv_results$warning[i] <- logistic$warning
  }
  # Pool the ten DIF estimates using Rubin's combining rules.
  m <- nrow(pv_results)
  logistic_estimate <- mean(pv_results$estimate)
  logistic_se <- sqrt(mean(pv_results$variance) + (1 + 1 / m) * var(pv_results$estimate))
  # Estimate DID from paired responses and report both methods on their own scales.
  did <- fit_did(dat$y_test, dat$y_anchor, dat$group)
  estimates <- c(logistic_estimate, unname(did["estimate"]))
  standard_errors <- c(logistic_se, unname(did["se"]))
  empirical <- data.frame(method = c("Logistic-regression DIF", "DID"),
    estimate = estimates, se = standard_errors,
    lower = estimates - Z * standard_errors, upper = estimates + Z * standard_errors,
    p = 2 * pnorm(abs(estimates / standard_errors), lower.tail = FALSE))

  # Use the PV mean difference as reference information, not as a fixed input.
  mu_reference <- mean(pv_results$mu)
  bounds <- c(0, EMPIRICAL_ETA * (4 * pnorm(abs(EMPIRICAL_MU) / 2) - 2),
              2 * EMPIRICAL_ETA)
  empirical_sensitivity <- data.frame(
    method = c("Unadjusted", "Normal", "Distribution-free"),
    eta_assumed = c(NA, EMPIRICAL_ETA, EMPIRICAL_ETA),
    mu_assumed = c(NA, EMPIRICAL_MU, NA), bound = bounds,
    lower = unname(did["estimate"] - Z * did["se"]) - bounds,
    upper = unname(did["estimate"] + Z * did["se"]) + bounds)
  reference_information <- data.frame(
    absolute_standardized_mean_difference = abs(mu_reference),
    eta_assumed = EMPIRICAL_ETA,
    absolute_mu_assumed = abs(EMPIRICAL_MU))
  write.csv(reference_information,
            file.path(folder, "sensitivity_inputs.csv"), row.names = FALSE)
  # Write the empirical estimates, sensitivity intervals, and supporting inputs.
  # Interval endpoints are calculated before rounding for presentation.
  report <- vapply(seq_len(nrow(empirical_sensitivity)), function(i) {
    x <- empirical_sensitivity[i, ]
    paste0(x$method, ": bound = ", fmt(x$bound),
           "; interval = [", fmt(x$lower), ", ", fmt(x$upper), "]")
  }, character(1))
  writeLines(report, file.path(folder, "sensitivity_report.txt"))
  write.csv(empirical, file.path(folder, "empirical_results.csv"), row.names = FALSE)
  write.csv(empirical_sensitivity, file.path(folder, "sensitivity.csv"), row.names = FALSE)
  write.csv(pv_results, file.path(folder, "pv_results.csv"), row.names = FALSE)
  writeLines(c(paste("Input:", basename(PIAAC_FILE)),
               paste("MD5:", unname(tools::md5sum(PIAAC_FILE))),
               paste("N:", nrow(dat)), paste("Reference:", sum(dat$group == 0)),
               paste("Focal:", sum(dat$group == 1)),
               "Unweighted; survey-design uncertainty not included.",
               "Sensitivity inputs treated as fixed; calibration uncertainty not included.",
               paste("Specified eta:", EMPIRICAL_ETA),
               paste("Specified mu:", EMPIRICAL_MU),
               paste("Reference absolute standardized mean difference:", abs(mu_reference))),
             file.path(folder, "analysis_notes.txt"))
  print(empirical, row.names = FALSE)
  print(empirical_sensitivity, row.names = FALSE)
}

run_empirical()

# 6. Figures 4 and 5 --------------------------------------------------------
# Shared settings give both figures the same text and line sizes at final width.

DEVICE_WIDTH <- 9.4
DEVICE_HEIGHT <- 8.85
FIGURE5_HEIGHT <- 10.5  # Eta rows and sample-size columns; coverage versus length.
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

# Center a legend with equal gaps between entries, using each label's actual width.
# Line-and-marker keys and marker-only keys use the same spacing rule.
draw_legend <- function(labels, pch, fills, colors = "black",
                        line_types = NULL, line_widths = NULL) {
  n <- length(labels)
  colors <- rep(colors, length.out = n)
  cex <- unname(TEXT_CEX["legend"])
  label_width <- strwidth(labels, cex = cex)
  key_width <- strwidth(if (is.null(line_types)) "M" else "MM", cex = cex)
  label_gap <- strwidth(" ", cex = cex)
  # Marker-only legends use wider gaps; line legends separate DID from Logistic.
  gap_multiplier <- if (is.null(line_types)) 2 else 1
  entry_gap <- strwidth("   ", cex = cex) * gap_multiplier
  entry_width <- key_width + label_gap + label_width
  extra_gap <- if (is.null(line_types)) 0 else entry_gap
  total_width <- sum(entry_width) + (n - 1) * entry_gap + extra_gap

  if (total_width > .98) {
    stop("Legend exceeds the figure width; review labels or text size.")
  }

  starts <- (1 - total_width) / 2 +
    c(0, head(cumsum(entry_width + entry_gap), -1))
  starts[-1] <- starts[-1] + extra_gap
  for (i in seq_len(n)) {
    if (!is.null(line_types)) {
      segments(starts[i], .5, starts[i] + key_width, .5,
               col = colors[i], lty = line_types[i], lwd = line_widths[i])
    }
    points(starts[i] + key_width / 2, .5, pch = pch[i],
           bg = fills[i], col = colors[i], cex = POINT_CEX,
           lwd = LINE_LWD["marker"])
    text(starts[i] + key_width + label_gap, .5, labels[i],
         adj = c(0, .5), cex = cex)
  }
}

# Figure 4: rejection rates under invalidity (A) and unreliability (B).
make_figure4 <- function() {
  folder <- file.path(OUTPUT_DIR, "study1")
  did <- read_summary("study1/did_summary.csv")
  logistic <- read_summary("study1/logistic_summary.csv")
  stopifnot(all(c("N", "tau", "rejection") %in% names(did)),
            all(c("N", "tau", "delta", "rho", "rejection") %in% names(logistic)))
  taus <- c(0, -.05, -.10)
  open_figure(file.path(folder, "figure4.pdf"))
  par(oma = c(.4, 4.8, .4, .3), family = FIGURE_FONT,
      ps = DEVICE_POINTSIZE, lwd = LINE_LWD["axis"])
  layout(
    rbind(c(1, 1, 1), c(2, 3, 4), c(5, 5, 5), c(0, 0, 0),
          c(6, 6, 6), c(7, 8, 9), c(10, 10, 10), c(11, 11, 11)),
    heights = c(.32, 2.3, .48, .10, .32, 2.3, .24, .48)
  )
  colors <- c("black", "grey35", "grey45", "grey55")
  line_types <- c(1, 2, 4, 5)
  symbols <- c(21, 21, 24, 22)
  widths <- unname(LINE_LWD[c("primary", "secondary", "secondary", "secondary")])
  for (block in 1:2) {
    blank_strip()
    text(0, .5, if (block == 1) "A" else "B", adj = c(0, .5), cex = TEXT_CEX["panel"], font = 2)
    text(
      .055, .5,
      if (block == 1) {
        "Matching-variable invalidity"
      } else {
        "Matching-variable unreliability"
      },
      adj = c(0, .5),
      cex = TEXT_CEX["title"]
    )
    for (n in c(500, 1000, 2000)) {
      par(mar = c(2.3, 2.8, 2.5, .8), mgp = c(1.7, .5, 0), tcl = -.2,
          las = 1, xaxs = "i", yaxs = "i", cex = 1)
      plot(NA, xlim = c(.005, -.105), ylim = c(-.025, 1.025),
           xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "l")
      abline(h = c(.25, .5, .75, 1), col = "grey90", lwd = LINE_LWD["grid"])
      abline(h = .05, col = "grey60", lty = 3, lwd = LINE_LWD["reference"])
      axis(1, at = taus, labels = c("0", "-0.05", "-0.10"), cex.axis = TEXT_CEX["tick"])
      axis(2, at = seq(0, 1, .2), labels = c("0", ".2", ".4", ".6", ".8", "1.0"), cex.axis = TEXT_CEX["tick"])
      if (block == 1L) {
        mtext(bquote(italic(N) == .(format(n, big.mark = ",", trim = TRUE))),
              side = 3, line = .7, at = mean(par("usr")[1:2]),
              adj = .5, cex = TEXT_CEX["title"])
      }
      z <- did[did$N == n, ]
      series <- list(z$rejection[match(taus, z$tau)])
      for (k in 1:3) {
        delta <- if (block == 1) c(0, .25, .5)[k] else 0
        rho <- if (block == 1) 1 else c(1, .8, .6)[k]
        z <- logistic[logistic$N == n & logistic$delta == delta & logistic$rho == rho, ]
        series[[k + 1]] <- z$rejection[match(taus, z$tau)]
      }
      stopifnot(all(is.finite(unlist(series))),
                all(unlist(series) >= 0 & unlist(series) <= 1))
      for (j in c(2, 3, 4, 1)) {
        lines(taus, series[[j]], col = colors[j], lty = line_types[j], lwd = widths[j])
      }
      for (j in c(2, 3, 4, 1)) {
        points(taus, series[[j]], col = colors[j], pch = symbols[j],
               bg = if (j == 1) "black" else "white", cex = POINT_CEX, lwd = LINE_LWD["marker"])
      }
    }
    if (block == 2L) {
      blank_strip()
      text(.5, .5, expression(tau), cex = TEXT_CEX["x"])
    }
    blank_strip()
    labels <- if (block == 1) {
      expression(DID, paste("Logistic: ", delta == 0),
                 paste("Logistic: ", delta == .25), paste("Logistic: ", delta == .50))
    } else {
      expression(DID, paste("Logistic: ", rho == 1.00),
                 paste("Logistic: ", rho == .80), paste("Logistic: ", rho == .60))
    }
    draw_legend(labels, symbols, c("black", rep("white", 3)),
                colors, line_types, widths)
  }
  mtext("Rejection rate", side = 2, outer = TRUE, line = 2.5,
        las = 0, at = .55, cex = TEXT_CEX["y"])
  dev.off()
}

make_figure4()

# Figure 5: joint coverage-length comparison.
# Rows are actual IRF gaps, matched by the specified eta; columns are N.
make_figure5 <- function() {
  folder <- file.path(OUTPUT_DIR, "study2")
  main <- read_summary("study2/main_summary.csv")
  eta_plot <- c(.05, .10, .15)
  n_plot <- c(500, 1000, 2000)
  method_plot <- c("Unadjusted", "Distribution-free", "Normal")
  method_labels <- c("Unadjusted", "Distribution-free bound", "Normal-distribution bound")
  method_symbols <- c(21, 24, 22)
  method_fills <- c("white", "white", "black")
  required_columns <- c("tau", "eta_true", "N", "method", "coverage", "mean_length")
  stopifnot(all(required_columns %in% names(main)))
  data <- main[abs(main$tau + .05) < 1e-8 & main$N %in% n_plot &
    main$method %in% method_plot & vapply(main$eta_true,
      function(x) any(abs(x - eta_plot) < 1e-8), logical(1)), required_columns]
  for (eta in eta_plot) for (n in n_plot) {
    z <- data[abs(data$eta_true - eta) < 1e-8 & data$N == n, ]
    stopifnot(nrow(z) == 3L, !anyDuplicated(z$method),
      setequal(z$method, method_plot), all(is.finite(z$coverage)),
      all(is.finite(z$mean_length)), all(z$coverage >= 0 & z$coverage <= 1),
      all(z$mean_length >= 0))
  }
  x_limits <- c(0, max(data$mean_length) * 1.08)
  y_limits <- c(max(0, floor(min(data$coverage) * 10) / 10 - .025), 1.035)
  x_ticks <- pretty(x_limits, n = 4)
  x_ticks <- x_ticks[x_ticks >= x_limits[1] & x_ticks <= x_limits[2]]
  y_ticks <- sort(unique(c(seq(ceiling(y_limits[1] * 10) / 10, 1, by = .1), .95)))
  x_labels <- sprintf("%.1f", x_ticks)
  x_labels[x_ticks == 0] <- "0"
  y_labels <- sub("^0\\.", ".", sprintf("%.2f", y_ticks))
  open_figure(file.path(folder, "figure5.pdf"), height = FIGURE5_HEIGHT)
  par(oma = c(.2, 4.5, .2, .3), family = FIGURE_FONT,
      ps = DEVICE_POINTSIZE, lwd = LINE_LWD["axis"])
  layout(rbind(c(1, 1, 1), c(2, 3, 4), c(5, 5, 5),
               c(6, 7, 8), c(9, 9, 9), c(10, 11, 12),
               c(13, 13, 13), c(14, 14, 14)),
         heights = c(.30, 2, .30, 2, .30, 2, .40, .55))
  for (i in seq_along(eta_plot)) {
    eta <- eta_plot[i]
    blank_strip()
    text(0, .5, LETTERS[i], adj = c(0, .5), font = 2, cex = TEXT_CEX["panel"])
    text(.055, .5, bquote(eta == .(sub("^0", "", sprintf("%.2f", eta)))),
         adj = c(0, .5), cex = TEXT_CEX["title"])
    for (j in seq_along(n_plot)) {
      z <- data[abs(data$eta_true - eta) < 1e-8 & data$N == n_plot[j], ]
      z <- z[match(method_plot, z$method), ]
      par(mar = c(2.2, 2.8, 2.5, .6), mgp = c(1.7, .5, 0),
          tcl = -.2, las = 1, xaxs = "i", yaxs = "i", cex = 1)
      plot(NA, xlim = x_limits, ylim = y_limits, xlab = "", ylab = "",
           xaxt = "n", yaxt = "n", bty = "l")
      if (i == 1L) {
        mtext(bquote(italic(N) == .(format(n_plot[j], big.mark = ",", trim = TRUE))),
              side = 3, line = .7, at = mean(par("usr")[1:2]),
              adj = .5, cex = TEXT_CEX["title"])
      }
      abline(h = setdiff(y_ticks, .95), col = "grey92", lwd = LINE_LWD["grid"])
      abline(h = .95, col = "grey50", lty = 3, lwd = LINE_LWD["reference"])
      axis(1, at = x_ticks, labels = FALSE)

      if (i == length(eta_plot)) {
        mtext(
          x_labels,
          side = 1,
          at = x_ticks,
          line = .6,
          cex = TEXT_CEX["tick"]
        )
      }
      axis(2, at = y_ticks, labels = if (j == 1L) y_labels else FALSE,
           cex.axis = TEXT_CEX["tick"])
      points(z$mean_length, z$coverage, pch = method_symbols, bg = method_fills,
             col = "black", cex = POINT_CEX, lwd = LINE_LWD["marker"])
    }
  }
  blank_strip()
  text(.5, .5, "Mean interval length", cex = TEXT_CEX["x"])
  blank_strip()
  draw_legend(method_labels, method_symbols, method_fills)
  mtext("Coverage", side = 2, outer = TRUE, line = 2.5,
        at = .56, las = 0, cex = TEXT_CEX["y"])
  dev.off()
}

make_figure5()

# 7. Supplementary tables ----------------------------------------------------

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
make_study1_tables <- function() {
  did_summary <- read_summary("study1/did_summary.csv")
  logistic_summary <- read_summary("study1/logistic_summary.csv")
  folder <- file.path(OUTPUT_DIR, "study1")
  d <- did_summary[order(did_summary$N, -did_summary$tau), ]
  tab <- data.frame(N = format(d$N, big.mark = ",", trim = TRUE),
    tau = fmt(d$tau, 2), bias = fmt(d$bias, 4), RMSE = fmt(d$rmse),
    SE_SD = fmt(d$se_sd_ratio), coverage = fmt(d$coverage, 4),
    rejection = fmt(d$rejection, 4))
  tab$N[duplicated(d$N)] <- ""
  pooled <- unique(d$n_total)
  stopifnot(length(pooled) == 1L)
  note <- paste0(
    "$N$ = sample size; $\\tau$ = true item-bias effect. ",
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
                    delta = fmt(keys$delta, 2))
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
    "$N$ = sample size; $\\tau$ = true item-bias effect; ",
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

make_study1_tables()

make_study2_tables <- function() {
  main <- read_summary("study2/main_summary.csv")
  joint <- read_summary("study2/eta_mu_summary.csv")
  stopifnot(all(abs(joint$tau - SENSITIVITY_TAU) < 1e-10),
            all(abs(joint$eta_true - SENSITIVITY_TRUE_ETA) < 1e-10),
            all(abs(joint$mu_true - MU) < 1e-10))
  folder <- file.path(OUTPUT_DIR, "study2")
  interval_note <- paste0(
    "Unadjusted denotes the nominal 95\\% DID confidence interval; ",
    "distribution-free and normal-distribution denote the intervals based ",
    "on the corresponding sensitivity bounds. Coverage is the proportion ",
    "of intervals containing $\\tau$; length is mean interval length.")
  repetitions <- unique(main$n_total)
  stopifnot(length(repetitions) == 1L)
  repetitions <- format(repetitions, big.mark = ",", trim = TRUE)

  # S3: additional effect sizes; specified eta equals the actual IRF gap.
  keys <- expand.grid(eta_true = ETA_TRUE_VALUES, N = N_VALUES, tau = c(0, -.10))
  tab <- data.frame(N = format(keys$N, big.mark = ",", trim = TRUE),
                    eta = fmt(keys$eta_true, 2))
  groups <- paste(keys$tau, keys$N)
  tab$N[duplicated(groups)] <- ""
  for (method in c("Unadjusted", "Distribution-free", "Normal")) {
    z <- main[main$method == method, ]
    index <- match(paste(keys$tau, keys$N, keys$eta_true),
                   paste(z$tau, z$N, z$eta_true))
    stopifnot(!anyNA(index))
    tab[[paste0("coverage", ncol(tab))]] <- fmt(z$coverage[index], 4)
    tab[[paste0("length", ncol(tab))]] <- fmt(z$mean_length[index])
  }
  note <- paste0(
    "$N$ = sample size; $\\eta$ = supremum of the absolute difference between ",
    "the test- and anchor-item IRFs under the reference-group condition; ",
    "$\\tau$ = true item-bias effect. ",
    interval_note, " Each condition is based on ", repetitions,
    " replications, with $\\widetilde{\\eta}=\\eta$ and ",
    "$\\widetilde{\\mu}=\\mu=-.5$. ",
    "The same simulated samples were used for all three intervals.")
  write_table(tab, c(
    " & & \\multicolumn{2}{c}{Unadjusted} & \\multicolumn{2}{c}{Distribution-free} & \\multicolumn{2}{c}{Normal-distribution} \\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
    "$N$ & $\\eta$ & Coverage & Length & Coverage & Length & Coverage & Length \\\\") ,
    "Interval Performance at Additional Effect Sizes in Study~2",
    "tab:sim2_other_effects", note, file.path(folder, "table_S3_other_effects.tex"),
    breaks = seq(4L, 16L, by = 3L),
    panel_titles = c("1" = "Panel A: $\\tau=0$", "10" = "Panel B: $\\tau=-.10$"))

  # S4: one baseline, three distribution-free settings, nine normal settings.
  # Sample sizes occupy columns; unused inputs are shown as dashes.
  keys <- rbind(
    data.frame(method = "Unadjusted", eta = ETA_SPECIFIED[1], mu = MU_SPECIFIED[1]),
    data.frame(method = "Distribution-free", eta = ETA_SPECIFIED, mu = MU_SPECIFIED[1]),
    data.frame(method = "Normal",
               expand.grid(mu = MU_SPECIFIED, eta = ETA_SPECIFIED,
                           KEEP.OUT.ATTRS = FALSE))[, c("method", "eta", "mu")])
  tab <- data.frame(eta = fmt(keys$eta, 3), mu = fmt(keys$mu, 2))
  tab$eta[keys$method == "Unadjusted"] <- "---"
  tab$mu[keys$method != "Normal"] <- "---"
  for (n in N_VALUES) {
    indices <- vapply(seq_len(nrow(keys)), function(i) {
      found <- which(joint$N == n & joint$method == keys$method[i] &
                     abs(joint$eta_assumed - keys$eta[i]) < 1e-10 &
                     abs(joint$mu_assumed - keys$mu[i]) < 1e-10)
      if (length(found) != 1L) stop("Missing or duplicate eta/mu setting in S4.")
      found
    }, integer(1))
    tab[[paste0("coverage_", n)]] <- fmt(joint$coverage[indices], 4)
    tab[[paste0("length_", n)]] <- fmt(joint$mean_length[indices])
  }
  repetitions <- unique(joint$n_total)
  stopifnot(length(repetitions) == 1L)
  note <- paste0(
    "$N$ = sample size; $\\widetilde{\\eta}$ and $\\widetilde{\\mu}$ = ",
    "specified sensitivity inputs; ",
    "$\\tau$ = true item-bias effect. ",
    "Panels A--C report results for the unadjusted nominal 95\\% DID confidence interval and ",
    "the intervals based on the distribution-free and normal-distribution bounds, respectively. ",
    "Coverage is the proportion of intervals containing $\\tau$; length is mean interval length. ",
    "A dash indicates a sensitivity parameter not used by the interval. ",
    "For each $N$, the same ", format(repetitions, big.mark = ",", trim = TRUE),
    " simulated samples were used across all three intervals and all 9 combinations, with $\\tau=-.05$, ",
    "the actual supremum of the absolute IRF difference $\\eta=.15$, ",
    "and the data-generating value $\\mu=-.5$.")
  write_table(tab, c(
    " & & \\multicolumn{2}{c}{$N=500$} & \\multicolumn{2}{c}{$N=1{,}000$} & \\multicolumn{2}{c}{$N=2{,}000$} \\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
    "$\\widetilde{\\eta}$ & $\\widetilde{\\mu}$ & Coverage & Length & Coverage & Length & Coverage & Length \\\\") ,
    "Sensitivity to the Specified Values of $\\widetilde{\\eta}$ and $\\widetilde{\\mu}$ in Study~2",
    "tab:sim2_eta_choice", note, file.path(folder, "table_S4_eta_mu_choice.tex"),
    breaks = c(2L, 5L, 8L, 11L),
    panel_titles = c("1" = "Panel A: Unadjusted", "2" = "Panel B: Distribution-free",
                     "5" = "Panel C: Normal-distribution"))
}

make_study2_tables()

# 8. Reproduction record ----------------------------------------------------
# Record the data file, analysis settings, and R session used for this run.
record_dir <- file.path(OUTPUT_DIR, "run_records")
dir.create(record_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(record_dir, "sessionInfo.txt"))
writeLines(capture.output(dput(list(
  piaac_file = basename(PIAAC_FILE), piaac_md5 = unname(tools::md5sum(PIAAC_FILE)),
  alpha = ALPHA, true_mu = MU, n_rep = N_REP, N = N_VALUES, tau = TAU_VALUES,
  eta_true = ETA_TRUE_VALUES, eta_specified = ETA_SPECIFIED, mu_specified = MU_SPECIFIED,
  empirical_eta = EMPIRICAL_ETA, empirical_mu = EMPIRICAL_MU,
  figure_font = FIGURE_FONT, final_figure_width_mm = FINAL_FIGURE_WIDTH_MM,
  intervals_clipped = FALSE, sensitivity_inputs_fixed = TRUE,
  study1_seed = 2026, study2_seed = 2027, seed_increment = 1009,
  study2_condition_ids = 10:36,
  run_time = format(Sys.time(), tz = "UTC")))),
  file.path(record_dir, "settings.txt"))
message("Finished. Results: ", normalizePath(OUTPUT_DIR))
