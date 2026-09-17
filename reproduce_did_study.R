# ==============================================================================
# Identifying Item Bias Without Conditioning:
# A Difference-in-Differences Approach
#
# Reproducibility script: Study 1, Study 2, and the empirical illustration
#
# CONTENTS
#   A. Common functions
#   B. Study 1: simulation
#   C. Study 1: tables and figures
#   D. Study 2: simulation
#   E. Study 2: tables and figures
#   F. Numerical checks for the simulation implementation
#   G. Empirical illustration: 2023 PIAAC Korean sample
#   H. Run the selected tasks
#
# REQUIREMENTS AND INPUT
#   - Simulations: R >= 3.6.0, using base R only.
#   - Empirical illustration: dplyr and sandwich, plus prgkorp2.csv.
#     Install packages once if needed: install.packages(c("dplyr", "sandwich"))
#   - Download the Korean public-use CSV from the OECD PIAAC 2nd Cycle Database.
#     Set EMPIRICAL_DATA_FILE below to its local path. Raw data are not included.
#   - Paths are relative to the current working directory unless absolute.
#
# HOW TO RUN
#   From this directory, use one of the following terminal commands:
#     Rscript reproduce_did_study.R all          # Full reproduction (default)
#     Rscript reproduce_did_study.R simulations  # Both simulation studies
#     Rscript reproduce_did_study.R empirical    # PIAAC illustration only
#     Rscript reproduce_did_study.R outputs      # Summarize saved simulations
#     Rscript reproduce_did_study.R checks       # Numerical checks only
#     Rscript reproduce_did_study.R smoke        # 20 replications per condition
#   In RStudio, set RUN_MODE below and source the entire file.
#   Use RUN_MODE = "custom" to choose individual tasks with the six switches.
#   The full simulation uses 5,000 replications per condition and may take
#   substantial time. The empirical analysis requires a local PIAAC CSV.
#
# OUTPUTS (under OUTPUT_ROOT)
#   study1/      Study 1 replications, Tables S1-S2, and Figure 4
#   study2/      Study 2 replications, Tables S3-S6, and Figure 5
#   empirical/   Empirical Table 1, PV summaries, sensitivity intervals,
#                metadata, and session information
#   Smoke mode writes to OUTPUT_ROOT/smoke/ and cannot replace full-run outputs.
#   Rerunning a selected task replaces its files in the selected output folder.
#
# ANALYTIC NOTES
#   Study 2 keeps the test-item IRF fixed across eta and uses
#     P_A(theta; eta) = eta + (1 - 2 * eta) * plogis(1.5 * theta).
#   Its saved model identifier prevents reuse of incompatible legacy results.
#   The empirical analysis is unweighted and treats eta and the estimated
#   standardized group difference as fixed sensitivity values.
#   Empirical calculations follow the supplied R 4.5.0 analysis using
#   dplyr 1.1.4 and sandwich 3.1-1; installed versions are saved for each run.
# ==============================================================================


# ==============================================================================
# EXECUTION SETTINGS
# ==============================================================================

# Choose: all, simulations, empirical, outputs, checks, smoke, or custom.
# A command-line mode overrides this value when run directly with Rscript.
RUN_MODE <- "all"

# These individual switches are used only in custom mode.
RUN_STUDY1_SIMULATION <- TRUE
RUN_STUDY1_OUTPUTS    <- TRUE
RUN_STUDY2_SIMULATION <- TRUE
RUN_STUDY2_OUTPUTS    <- TRUE
RUN_CHECKS           <- TRUE
RUN_EMPIRICAL        <- TRUE

N_REP       <- 5000L
STUDY1_SEED <- 2026L
STUDY2_SEED <- 2027L
OUTPUT_ROOT <- "results"

# Empirical analysis: local input and fixed IRF-nonequivalence bound.
EMPIRICAL_DATA_FILE <- "prgkorp2.csv"
EMPIRICAL_ETA       <- .033

ALPHA <- .05
STUDY2_MODEL_ID <- "fixed_test_smooth_anchor_v2"

# Apply the selected execution mode. No packages or data are needed for
# simulations, saved simulation outputs, or the base-R numerical checks.
if (sys.nframe() == 0L) {
  run_arguments <- commandArgs(trailingOnly = TRUE)
  if (length(run_arguments) > 1L) {
    stop("Supply at most one mode, for example: Rscript reproduce_did_study.R checks")
  }
  if (length(run_arguments) == 1L) {
    RUN_MODE <- run_arguments[1L]
  }
}

allowed_modes <- c(
  "all", "simulations", "empirical", "outputs", "checks", "smoke", "custom"
)
if (length(RUN_MODE) != 1L || is.na(RUN_MODE) ||
    !RUN_MODE %in% allowed_modes) {
  stop("RUN_MODE must be one of: ", paste(allowed_modes, collapse = ", "))
}

if (RUN_MODE != "custom") {
  RUN_STUDY1_SIMULATION <- RUN_MODE %in% c("all", "simulations", "smoke")
  RUN_STUDY2_SIMULATION <- RUN_STUDY1_SIMULATION
  RUN_STUDY1_OUTPUTS <- RUN_MODE %in% c("all", "simulations", "outputs", "smoke")
  RUN_STUDY2_OUTPUTS <- RUN_STUDY1_OUTPUTS
  RUN_CHECKS <- RUN_MODE %in% c("all", "simulations", "checks", "smoke")
  RUN_EMPIRICAL <- RUN_MODE %in% c("all", "empirical")
}

if (RUN_MODE == "smoke") {
  N_REP <- 20L
  OUTPUT_ROOT <- file.path(OUTPUT_ROOT, "smoke")
  message("Smoke mode checks execution only; its results are not manuscript results.")
}


# ==============================================================================
# A. COMMON FUNCTIONS
# ==============================================================================

if (getRversion() < "3.6.0") {
  stop("This script requires R >= 3.6.0.")
}


write_csv <- function(x, filename) {
  write.csv(x, file = filename, row.names = FALSE, na = "")
  invisible(filename)
}


prepare_directory <- function(directory) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(directory)) {
    stop("Could not create directory: ", directory)
  }
  invisible(directory)
}


set_condition_seed <- function(base_seed, condition_id) {
  RNGkind(
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )

  seed <- as.double(base_seed) + 1009 * as.double(condition_id)

  if (!is.finite(seed) ||
      seed < 1 ||
      seed > .Machine$integer.max) {
    stop("Invalid condition-specific seed.")
  }

  set.seed(as.integer(seed))
  invisible(seed)
}


save_run_information <- function(metadata, directory) {
  writeLines(
    capture.output(dput(metadata)),
    file.path(directory, "metadata.txt")
  )

  writeLines(
    capture.output(sessionInfo()),
    file.path(directory, "sessionInfo.txt")
  )

  invisible(NULL)
}


make_metadata <- function(study, replications, base_seed, alpha) {
  list(
    study = study,
    replications = as.integer(replications),
    base_seed = as.integer(base_seed),
    condition_seed_rule = "base_seed + 1009 * condition_id",
    alpha = alpha,
    rng_kind = "Mersenne-Twister",
    normal_kind = "Inversion",
    sample_kind = "Rejection",
    group_probability = .5,
    mu = -.5,
    latent_sd = 1,
    anchor_discrimination = 1.5,
    anchor_difficulty = 0,
    responses_conditionally_independent = TRUE,
    sensitivity_inputs_fixed = TRUE,
    confidence_intervals_clipped = FALSE,
    R_version = R.version.string,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}


integrate_real <- function(fun) {
  integrate(
    fun,
    lower = -Inf,
    upper = Inf,
    subdivisions = 1000L,
    rel.tol = 1e-10,
    abs.tol = 1e-12,
    stop.on.error = TRUE
  )$value
}


anchor_irf <- function(theta) {
  plogis(1.5 * theta)
}


# Study 1 retains the baseline anchor; Study 2 perturbs only the anchor.
comparison_anchor_irf <- function(theta, violation, parameter) {
  if (identical(violation, "Smooth")) {
    if (length(parameter) != 1L || !is.finite(parameter) ||
        parameter < 0 || parameter >= .5) {
      stop("Smooth nonequivalence requires 0 <= eta < .5.")
    }
    return(parameter + (1 - 2 * parameter) * anchor_irf(theta))
  }
  anchor_irf(theta)
}


reference_lp <- function(theta, violation, parameter) {
  if (identical(violation, "Smooth")) {
    if (length(parameter) != 1L ||
        !is.finite(parameter) ||
        parameter < 0 ||
        parameter >= .5) {
      stop("Smooth nonequivalence requires 0 <= eta < .5.")
    }

    # P_T(0, theta) stays fixed for every eta.
    return(1.5 * theta)
  }

  if (identical(violation, "Difficulty")) {
    return(1.5 * (theta - parameter))
  }

  if (identical(violation, "Discrimination")) {
    return(parameter * theta)
  }

  stop("Unknown violation type: ", violation)
}


reference_irf <- function(theta, violation, parameter) {
  plogis(reference_lp(theta, violation, parameter))
}


max_difference_difficulty <- function(delta_b) {
  if (length(delta_b) != 1L || !is.finite(delta_b)) {
    stop("delta_b must be one finite number.")
  }

  tanh(1.5 * abs(delta_b) / 4)
}


max_difference_discrimination <- function(a_test) {
  a_anchor <- 1.5

  if (length(a_test) != 1L ||
      !is.finite(a_test) ||
      a_test <= 0 ||
      a_test > a_anchor) {
    stop("a_test must satisfy 0 < a_test <= 1.5.")
  }

  if (a_test == a_anchor) {
    return(0)
  }

  # At the positive maximizer:
  # a_test * logistic'(a_test * theta)
  #   = a_anchor * logistic'(a_anchor * theta).
  # Log derivatives avoid numerical underflow.
  derivative_log_difference <- function(theta) {
    log(a_test) +
      dlogis(a_test * theta, log = TRUE) -
      log(a_anchor) -
      dlogis(a_anchor * theta, log = TRUE)
  }

  upper <- 1

  while (derivative_log_difference(upper) <= 0) {
    upper <- 2 * upper

    if (!is.finite(upper) || upper > 1e8) {
      stop("Could not bracket the discrimination-gap maximizer.")
    }
  }

  theta_max <- uniroot(
    derivative_log_difference,
    interval = c(0, upper),
    tol = 1e-12
  )$root

  plogis(a_anchor * theta_max) - plogis(a_test * theta_max)
}


calibrate_parameter <- function(eta, violation) {
  if (length(eta) != 1L ||
      !is.finite(eta) ||
      eta < 0) {
    stop("eta must be one finite nonnegative number.")
  }

  if (identical(violation, "Smooth")) {
    if (eta >= .5) {
      stop("Smooth nonequivalence requires eta < .5.")
    }

    return(eta)
  }

  if (identical(violation, "Difficulty")) {
    if (eta >= 1) {
      stop("Difficulty calibration requires eta < 1.")
    }

    return(4 * atanh(eta) / 1.5)
  }

  if (identical(violation, "Discrimination")) {
    if (eta >= .5) {
      stop("This discrimination calibration requires eta < .5.")
    }

    if (eta == 0) {
      return(1.5)
    }

    lower <- 1e-6

    while (max_difference_discrimination(lower) < eta) {
      lower <- lower / 10

      if (lower < 1e-14) {
        stop("Could not bracket the discrimination calibration.")
      }
    }

    return(
      uniroot(
        function(a_test) {
          max_difference_discrimination(a_test) - eta
        },
        interval = c(lower, 1.5),
        tol = 1e-12
      )$root
    )
  }

  stop("Unknown violation type: ", violation)
}


tau_from_gamma <- function(
    gamma,
    violation,
    parameter,
    mu = -.5) {

  integrate_real(function(theta) {
    lp <- reference_lp(theta, violation, parameter)

    (plogis(lp + gamma) - plogis(lp)) *
      dnorm(theta, mean = mu, sd = 1)
  })
}


calibrate_gamma <- function(
    tau,
    violation,
    parameter,
    mu = -.5) {

  if (length(tau) != 1L || !is.finite(tau)) {
    stop("tau must be one finite number.")
  }

  if (tau == 0) {
    return(0)
  }

  mean_reference <- integrate_real(function(theta) {
    reference_irf(theta, violation, parameter) *
      dnorm(theta, mean = mu, sd = 1)
  })

  if (tau <= -mean_reference ||
      tau >= 1 - mean_reference) {
    stop("Requested tau is outside the attainable open interval.")
  }

  objective <- function(gamma) {
    tau_from_gamma(gamma, violation, parameter, mu) - tau
  }

  lower <- -1
  upper <- 1

  while (objective(lower) > 0) {
    lower <- 2 * lower

    if (abs(lower) > 1e6) {
      stop("Could not bracket gamma below.")
    }
  }

  while (objective(upper) < 0) {
    upper <- 2 * upper

    if (abs(upper) > 1e6) {
      stop("Could not bracket gamma above.")
    }
  }

  uniroot(
    objective,
    interval = c(lower, upper),
    tol = 1e-12
  )$root
}


identification_error <- function(
    violation,
    parameter,
    mu = -.5) {

  if (mu == 0) {
    return(0)
  }

  integrate_real(function(theta) {
    gap <- reference_irf(theta, violation, parameter) -
      comparison_anchor_irf(theta, violation, parameter)

    density_difference <- dnorm(theta, mean = mu, sd = 1) -
      dnorm(theta, mean = 0, sd = 1)

    gap * density_difference
  })
}


population_did <- function(
    gamma,
    violation,
    parameter,
    mu = -.5) {

  focal_difference <- integrate_real(function(theta) {
    (
      plogis(reference_lp(theta, violation, parameter) + gamma) -
        comparison_anchor_irf(theta, violation, parameter)
    ) * dnorm(theta, mean = mu, sd = 1)
  })

  reference_difference <- integrate_real(function(theta) {
    (
      reference_irf(theta, violation, parameter) -
        comparison_anchor_irf(theta, violation, parameter)
    ) * dnorm(theta, mean = 0, sd = 1)
  })

  focal_difference - reference_difference
}


normal_bound <- function(eta, mu) {
  if (any(!is.finite(eta)) ||
      any(eta < 0) ||
      any(!is.finite(mu))) {
    stop("Invalid inputs to normal_bound().")
  }

  eta * (4 * pnorm(abs(mu) / 2) - 2)
}


distribution_free_bound <- function(eta) {
  if (any(!is.finite(eta)) || any(eta < 0)) {
    stop("Invalid eta in distribution_free_bound().")
  }

  2 * eta
}


fit_did_hc3 <- function(y_test, y_anchor, group) {
  if (length(y_test) != length(y_anchor) ||
      length(y_test) != length(group)) {
    stop("Response and group vectors must have equal lengths.")
  }

  if (anyNA(y_test) ||
      anyNA(y_anchor) ||
      anyNA(group) ||
      any(!is.finite(y_test)) ||
      any(!is.finite(y_anchor)) ||
      !all(group %in% c(0, 1))) {
    stop("Invalid data in fit_did_hc3().")
  }

  difference <- y_test - y_anchor
  d0 <- difference[group == 0]
  d1 <- difference[group == 1]

  n0 <- length(d0)
  n1 <- length(d1)

  if (n0 < 2L || n1 < 2L) {
    stop("HC3 DID estimation requires at least two observations per group.")
  }

  estimate <- mean(d1) - mean(d0)

  # Exact HC3 variance for the intercept-plus-binary-group regression.
  variance <- var(d1) / (n1 - 1) +
    var(d0) / (n0 - 1)

  c(
    estimate = estimate,
    se = sqrt(variance),
    n0 = n0,
    n1 = n1
  )
}


fit_logistic_dif <- function(y_test, matching, group) {
  warning_messages <- character(0)
  error_message <- ""

  fit <- tryCatch(
    withCallingHandlers(
      glm(
        y_test ~ matching + group,
        family = binomial(link = "logit"),
        control = glm.control(maxit = 50L)
      ),
      warning = function(w) {
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )

  invalid_result <- function(reason) {
    list(
      estimate = NA_real_,
      se = NA_real_,
      valid = FALSE,
      failed = TRUE,
      warning_count = length(warning_messages),
      warning = length(warning_messages) > 0L,
      reason = reason
    )
  }

  if (is.null(fit)) {
    return(invalid_result(
      if (nzchar(error_message)) error_message else "glm_error"
    ))
  }

  if (!isTRUE(fit$converged)) {
    return(invalid_result("nonconvergence"))
  }

  if (isTRUE(fit$boundary)) {
    return(invalid_result("boundary_fit"))
  }

  coefficient_table <- tryCatch(
    coef(summary(fit)),
    error = function(e) NULL
  )

  if (is.null(coefficient_table) ||
      !("group" %in% rownames(coefficient_table))) {
    return(invalid_result("missing_group_coefficient"))
  }

  estimate <- unname(coefficient_table["group", "Estimate"])
  se <- unname(coefficient_table["group", "Std. Error"])

  if (!is.finite(estimate) ||
      !is.finite(se) ||
      se <= 0) {
    return(invalid_result("invalid_estimate_or_se"))
  }

  list(
    estimate = estimate,
    se = se,
    valid = TRUE,
    failed = FALSE,
    warning_count = length(warning_messages),
    warning = length(warning_messages) > 0L,
    reason = ""
  )
}


binomial_mcse <- function(probability, n) {
  if (!is.finite(n) || n <= 0) {
    return(rep(NA_real_, length(probability)))
  }

  sqrt(probability * (1 - probability) / n)
}


performance <- function(
    estimates,
    standard_errors,
    truth,
    alpha = .05) {

  if (length(estimates) != length(standard_errors)) {
    stop("Estimate and SE vectors must have equal lengths.")
  }

  if (length(truth) != 1L ||
      !is.finite(truth) ||
      length(alpha) != 1L ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1) {
    stop("Invalid truth or alpha in performance().")
  }

  valid <- is.finite(estimates) &
    is.finite(standard_errors) &
    standard_errors >= 0

  n_total <- length(estimates)
  n_valid <- sum(valid)

  if (n_valid < 2L) {
    stop("At least two valid replications are required.")
  }

  estimates <- estimates[valid]
  standard_errors <- standard_errors[valid]

  critical <- qnorm(1 - alpha / 2)
  lower <- estimates - critical * standard_errors
  upper <- estimates + critical * standard_errors

  empirical_sd <- sd(estimates)
  mean_se <- mean(standard_errors)

  coverage <- mean(lower <= truth & upper >= truth)
  rejection <- mean(lower > 0 | upper < 0)
  negative <- mean(upper < 0)
  positive <- mean(lower > 0)

  data.frame(
    n_total = n_total,
    n_valid = n_valid,
    n_invalid = n_total - n_valid,
    mean_estimate = mean(estimates),
    bias = mean(estimates) - truth,
    bias_mcse = empirical_sd / sqrt(n_valid),
    empirical_sd = empirical_sd,
    mean_se = mean_se,
    se_sd_ratio = if (empirical_sd > 0) {
      mean_se / empirical_sd
    } else {
      NA_real_
    },
    rmse = sqrt(mean((estimates - truth)^2)),
    coverage = coverage,
    coverage_mcse = binomial_mcse(coverage, n_valid),
    mean_width = mean(upper - lower),
    rejection = rejection,
    rejection_mcse = binomial_mcse(rejection, n_valid),
    excludes_zero_negative = negative,
    excludes_zero_positive = positive,
    stringsAsFactors = FALSE
  )
}


sensitivity_summary <- function(
    estimates,
    standard_errors,
    tau,
    b_normal,
    b_df,
    alpha = .05) {

  if (length(estimates) != length(standard_errors)) {
    stop("Estimate and SE vectors must have equal lengths.")
  }

  if (length(tau) != 1L ||
      !is.finite(tau) ||
      length(b_normal) != 1L ||
      length(b_df) != 1L ||
      !is.finite(b_normal) ||
      !is.finite(b_df) ||
      b_normal < 0 ||
      b_df < b_normal - 1e-12 ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1) {
    stop("Invalid inputs to sensitivity_summary().")
  }

  valid <- is.finite(estimates) &
    is.finite(standard_errors) &
    standard_errors >= 0

  n_total <- length(estimates)
  n_valid <- sum(valid)

  if (n_valid < 2L) {
    stop("At least two valid replications are required.")
  }

  estimates <- estimates[valid]
  standard_errors <- standard_errors[valid]

  critical <- qnorm(1 - alpha / 2)

  unadjusted_lower <- estimates - critical * standard_errors
  unadjusted_upper <- estimates + critical * standard_errors

  methods <- c(
    "Unadjusted",
    "Normal",
    "Distribution_free"
  )

  bounds <- c(0, b_normal, b_df)

  answer <- lapply(seq_along(methods), function(j) {
    # No clipping to [-1, 1].
    lower <- unadjusted_lower - bounds[j]
    upper <- unadjusted_upper + bounds[j]

    covered <- lower <= tau & upper >= tau
    excluded_negative <- upper < 0
    excluded_positive <- lower > 0
    excluded <- excluded_negative | excluded_positive

    coverage <- mean(covered)
    exclusion <- mean(excluded)
    negative <- mean(excluded_negative)
    positive <- mean(excluded_positive)

    widths <- upper - lower

    data.frame(
      method = methods[j],
      bound = bounds[j],
      n_total = n_total,
      n_valid = n_valid,
      n_invalid = n_total - n_valid,
      coverage_tau = coverage,
      coverage_mcse = binomial_mcse(coverage, n_valid),
      mean_width = mean(widths),
      mean_width_mcse = sd(widths) / sqrt(n_valid),
      excludes_zero = exclusion,
      excludes_zero_mcse = binomial_mcse(exclusion, n_valid),
      excludes_zero_negative = negative,
      excludes_zero_negative_mcse = binomial_mcse(negative, n_valid),
      excludes_zero_positive = positive,
      excludes_zero_positive_mcse = binomial_mcse(positive, n_valid),
      stringsAsFactors = FALSE
    )
  })

  answer <- do.call(rbind, answer)
  rownames(answer) <- NULL
  answer
}


# ==============================================================================
# B. STUDY 1 -- SIMULATION
# ==============================================================================

simulate_study1 <- function(
    replications,
    base_seed,
    directory) {

  prepare_directory(directory)

  if (file.exists(file.path(directory, "simulation.rds"))) {
    message("Study 1 simulation is enabled: existing saved results will be overwritten.")
  }

  mu <- -.5

  design <- expand.grid(
    N = c(500L, 1000L, 2000L),
    tau = c(0, -.05, -.10),
    delta = c(0, .25, .50),
    rho_M = c(1, .8, .6),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  design$condition_id <- seq_len(nrow(design))
  design$mu <- mu

  calibration <- data.frame(
    tau = c(0, -.05, -.10),
    stringsAsFactors = FALSE
  )

  calibration$gamma <- vapply(
    calibration$tau,
    function(tau) {
      calibrate_gamma(
        tau = tau,
        violation = "Difficulty",
        parameter = 0,
        mu = mu
      )
    },
    numeric(1)
  )

  calibration$tau_achieved <- vapply(
    seq_len(nrow(calibration)),
    function(i) {
      tau_from_gamma(
        gamma = calibration$gamma[i],
        violation = "Difficulty",
        parameter = 0,
        mu = mu
      )
    },
    numeric(1)
  )

  design$gamma <- calibration$gamma[
    match(design$tau, calibration$tau)
  ]

  design$condition_seed <- as.double(base_seed) +
    1009 * design$condition_id

  design <- design[
    ,
    c(
      "condition_id", "N", "tau", "delta", "rho_M",
      "mu", "gamma", "condition_seed"
    )
  ]

  stopifnot(
    nrow(design) == 81L,
    max(abs(calibration$tau_achieved - calibration$tau)) < 1e-8
  )

  raw_list <- vector("list", nrow(design))

  for (i in seq_len(nrow(design))) {
    d <- design[i, ]
    set_condition_seed(base_seed, d$condition_id)

    did <- numeric(replications)
    did_se <- numeric(replications)
    n0 <- integer(replications)
    n1 <- integer(replications)

    logistic <- rep(NA_real_, replications)
    logistic_se <- rep(NA_real_, replications)
    logistic_valid <- logical(replications)
    logistic_failed <- logical(replications)
    logistic_warning <- logical(replications)
    logistic_warning_count <- integer(replications)
    logistic_failure_reason <- character(replications)

    matching_error_sd <- sqrt((1 - d$rho_M) / d$rho_M)

    for (r in seq_len(replications)) {
      group <- rbinom(d$N, size = 1L, prob = .5)
      theta <- rnorm(d$N, mean = mu * group, sd = 1)

      y_anchor <- rbinom(
        d$N,
        size = 1L,
        prob = anchor_irf(theta)
      )

      y_test <- rbinom(
        d$N,
        size = 1L,
        prob = plogis(1.5 * theta + d$gamma * group)
      )

      matching_error <- if (matching_error_sd == 0) {
        numeric(d$N)
      } else {
        rnorm(d$N, mean = 0, sd = matching_error_sd)
      }

      matching <- theta + d$delta * group + matching_error

      did_fit <- fit_did_hc3(y_test, y_anchor, group)
      logistic_fit <- fit_logistic_dif(y_test, matching, group)

      did[r] <- did_fit["estimate"]
      did_se[r] <- did_fit["se"]
      n0[r] <- did_fit["n0"]
      n1[r] <- did_fit["n1"]

      logistic[r] <- logistic_fit$estimate
      logistic_se[r] <- logistic_fit$se
      logistic_valid[r] <- logistic_fit$valid
      logistic_failed[r] <- logistic_fit$failed
      logistic_warning[r] <- logistic_fit$warning
      logistic_warning_count[r] <- logistic_fit$warning_count
      logistic_failure_reason[r] <- logistic_fit$reason
    }

    raw_list[[i]] <- data.frame(
      condition_id = rep(d$condition_id, replications),
      replication = seq_len(replications),
      did = did,
      did_se = did_se,
      n0 = n0,
      n1 = n1,
      logistic = logistic,
      logistic_se = logistic_se,
      logistic_valid = logistic_valid,
      logistic_failed = logistic_failed,
      logistic_warning = logistic_warning,
      logistic_warning_count = logistic_warning_count,
      logistic_failure_reason = logistic_failure_reason,
      stringsAsFactors = FALSE
    )

    cat(
      sprintf(
        "Study 1: condition %d/%d completed (N=%d, tau=%.2f, delta=%.2f, rho_M=%.2f).\n",
        i, nrow(design), d$N, d$tau, d$delta, d$rho_M
      )
    )

    flush.console()
  }

  raw <- do.call(rbind, raw_list)
  rownames(raw) <- NULL

  metadata <- make_metadata(
    study = "Study 1",
    replications = replications,
    base_seed = base_seed,
    alpha = ALPHA
  )

  metadata$n_conditions <- nrow(design)
  metadata$matching_error_variance <- "(1 - rho_M) / rho_M"
  metadata$logistic_model <- "Y_T ~ M + G"
  metadata$logistic_inference <- "Model-based Wald confidence interval"

  simulation <- list(
    raw = raw,
    design = design,
    calibration = calibration,
    metadata = metadata
  )

  saveRDS(
    simulation,
    file.path(directory, "simulation.rds")
  )

  write_csv(
    design,
    file.path(directory, "study1_design.csv")
  )

  write_csv(
    calibration,
    file.path(directory, "study1_calibration.csv")
  )

  save_run_information(metadata, directory)

  cat("Study 1 simulation saved.\n")
  invisible(simulation)
}


# ==============================================================================
# C. STUDY 1 -- OUTPUTS ONLY
# ==============================================================================

outputs_study1 <- function(directory) {
  input_file <- file.path(directory, "simulation.rds")

  if (!file.exists(input_file)) {
    stop("No saved Study 1 simulation. Set RUN_STUDY1_SIMULATION <- TRUE first.")
  }

  simulation <- readRDS(input_file)

  required_raw <- c(
    "condition_id", "did", "did_se", "logistic", "logistic_se",
    "logistic_valid", "logistic_failed", "logistic_warning",
    "logistic_warning_count"
  )

  required_design <- c(
    "condition_id", "N", "tau", "delta", "rho_M", "gamma"
  )

  if (!all(required_raw %in% names(simulation$raw)) ||
      !all(required_design %in% names(simulation$design))) {
    stop(
      "Saved Study 1 data have a different schema. ",
      "Study 1 outputs were not changed; keep RUN_STUDY1_OUTPUTS <- FALSE."
    )
  }

  results <- merge(
    simulation$raw,
    simulation$design,
    by = "condition_id",
    sort = FALSE
  )

  alpha <- simulation$metadata$alpha
  groups <- split(results, results$condition_id)

  keys <- c(
    "condition_id", "N", "tau", "delta", "rho_M", "gamma"
  )

  cat(
    "Study 1 outputs use saved replications per condition: ",
    simulation$metadata$replications,
    "\n",
    sep = ""
  )

  did_summary <- do.call(rbind, lapply(groups, function(d) {
    cbind(
      d[1, keys],
      performance(
        estimates = d$did,
        standard_errors = d$did_se,
        truth = d$tau[1],
        alpha = alpha
      )
    )
  }))

  rownames(did_summary) <- NULL

  logistic_summary <- do.call(rbind, lapply(groups, function(d) {
    valid <- d$logistic_valid &
      is.finite(d$logistic) &
      is.finite(d$logistic_se) &
      d$logistic_se > 0

    n_valid <- sum(valid)
    n_total <- nrow(d)

    if (n_valid > 0L) {
      critical <- qnorm(1 - alpha / 2)
      lower <- d$logistic[valid] -
        critical * d$logistic_se[valid]
      upper <- d$logistic[valid] +
        critical * d$logistic_se[valid]

      rejection <- mean(lower > 0 | upper < 0)
      negative <- mean(upper < 0)
      positive <- mean(lower > 0)

      mean_coefficient <- mean(d$logistic[valid])
      coefficient_sd <- if (n_valid > 1L) {
        sd(d$logistic[valid])
      } else {
        NA_real_
      }

      mean_se <- mean(d$logistic_se[valid])
    } else {
      rejection <- NA_real_
      negative <- NA_real_
      positive <- NA_real_
      mean_coefficient <- NA_real_
      coefficient_sd <- NA_real_
      mean_se <- NA_real_
    }

    cbind(
      d[1, keys],
      data.frame(
        n_total = n_total,
        n_valid = n_valid,
        n_failed = n_total - n_valid,
        n_reported_failures = sum(d$logistic_failed),
        n_fits_with_warnings = sum(d$logistic_warning),
        n_warnings = sum(d$logistic_warning_count),
        mean_coefficient = mean_coefficient,
        coefficient_sd = coefficient_sd,
        mean_se = mean_se,
        rejection = rejection,
        rejection_mcse = binomial_mcse(rejection, n_valid),
        excludes_zero_negative = negative,
        excludes_zero_positive = positive,
        stringsAsFactors = FALSE
      )
    )
  }))

  rownames(logistic_summary) <- NULL

  # Pool raw DID replications, rather than averaging condition-level
  # nonlinear statistics such as empirical SD or RMSE.
  pooled_keys <- unique(results[, c("N", "tau")])
  pooled_keys <- pooled_keys[
    order(pooled_keys$N, -pooled_keys$tau),
    ,
    drop = FALSE
  ]

  did_pooled <- do.call(
    rbind,
    lapply(seq_len(nrow(pooled_keys)), function(i) {
      current_N <- pooled_keys$N[i]
      current_tau <- pooled_keys$tau[i]

      d <- results[
        results$N == current_N & results$tau == current_tau,
        ,
        drop = FALSE
      ]

      cbind(
        pooled_keys[i, , drop = FALSE],
        data.frame(
          n_conditions_pooled = length(unique(d$condition_id))
        ),
        performance(
          estimates = d$did,
          standard_errors = d$did_se,
          truth = current_tau,
          alpha = alpha
        )
      )
    })
  )

  rownames(did_pooled) <- NULL

  figure_data <- logistic_summary[
    logistic_summary$tau == 0,
    ,
    drop = FALSE
  ]

  # Assign pooled null DID values explicitly by sample size.
  did_null <- did_pooled[
    did_pooled$tau == 0,
    ,
    drop = FALSE
  ]

  figure_data$did_rejection <- did_null$rejection[
    match(figure_data$N, did_null$N)
  ]

  write_csv(
    did_summary,
    file.path(directory, "study1_did_summary.csv")
  )

  write_csv(
    did_pooled,
    file.path(directory, "study1_did_pooled_summary.csv")
  )

  write_csv(
    logistic_summary,
    file.path(directory, "study1_logistic_summary.csv")
  )

  write_csv(
    did_pooled,
    file.path(directory, "table_study1_did.csv")
  )

  write_csv(
    logistic_summary,
    file.path(directory, "table_study1_logistic.csv")
  )

  write_csv(
    figure_data,
    file.path(directory, "figure4_data.csv")
  )

  plot_figure <- function() {
    pdf(
      file.path(directory, "figure4.pdf"),
      width = 8.5,
      height = 4.3,
      family = "serif",
      useDingbats = FALSE
    )

    on.exit(dev.off(), add = TRUE)

    layout(
      matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE),
      heights = c(3.2, .72)
    )

    par(
      oma = c(0, 2.4, .2, .2),
      family = "serif"
    )

    rho_values <- c(1, .8, .6)
    rho_lty <- c(1, 2, 4)
    rho_pch <- c(21, 22, 24)

    for (current_N in c(500, 1000, 2000)) {
      par(
        mar = c(2.5, 2.8, 2.1, .8),
        mgp = c(1.6, .50, 0),
        tcl = -.23,
        las = 1,
        cex = 1,
        cex.axis = 1.05,
        cex.main = 1.12,
        xaxs = "i",
        yaxs = "i"
      )

      plot(
        NA,
        xlim = c(-.015, .515),
        ylim = c(0, 1.025),
        xlab = "",
        ylab = "",
        xaxt = "n",
        yaxt = "n",
        bty = "l",
        xaxs = "i",
        yaxs = "i"
      )

      title(
        main = bquote(
          italic(N) == .(
            format(
              current_N,
              big.mark = ",",
              trim = TRUE
            )
          )
        ),
        line = .65
      )

      axis(
        1,
        at = c(0, .25, .50),
        labels = c("0", ".25", ".50")
      )

      axis(
        2,
        at = seq(0, 1, by = .2),
        labels = sub("^0", "", sprintf("%.1f", seq(0, 1, by = .2)))
      )

      abline(
        h = alpha,
        col = "grey55",
        lty = 3,
        lwd = .8
      )

      for (j in seq_along(rho_values)) {
        d <- figure_data[
          figure_data$N == current_N &
            figure_data$rho_M == rho_values[j],
          ,
          drop = FALSE
        ]

        d <- d[order(d$delta), , drop = FALSE]

        lines(
          d$delta,
          d$rejection,
          lty = rho_lty[j],
          lwd = 1.30
        )
      }

      did_rate <- did_null$rejection[did_null$N == current_N]

      lines(
        c(0, .25, .50),
        rep(did_rate, 3L),
        lty = 3,
        lwd = 1.30
      )

      for (j in seq_along(rho_values)) {
        d <- figure_data[
          figure_data$N == current_N &
            figure_data$rho_M == rho_values[j],
          ,
          drop = FALSE
        ]

        d <- d[order(d$delta), , drop = FALSE]

        points(
          d$delta,
          d$rejection,
          pch = rho_pch[j],
          bg = "white",
          cex = .92,
          lwd = 1.05
        )
      }

      points(
        c(0, .25, .50),
        rep(did_rate, 3L),
        pch = 23,
        bg = "white",
        cex = .92,
        lwd = 1.05
      )
    }

    par(mar = c(0, 0, 0, 0))
    plot.new()

    text(
      .5,
      .84,
      labels = expression(
        paste("Matching-variable invalidity (", delta, ")")
      ),
      cex = 1.15
    )

    legend(
      x = .40,
      y = .37,
      xjust = .5,
      yjust = .5,
      legend = c(
        "Logistic: reliability 1.0",
        "Logistic: reliability .8",
        "Logistic: reliability .6"
      ),
      lty = rho_lty,
      pch = rho_pch,
      pt.bg = "white",
      pt.cex = .95,
      pt.lwd = 1.05,
      col = "black",
      horiz = TRUE,
      bty = "n",
      cex = .91,
      lwd = 1.30,
      seg.len = 2.0,
      x.intersp = .65
    )

    legend(
      x = .85,
      y = .37,
      xjust = .5,
      yjust = .5,
      legend = "DID",
      lty = 3,
      pch = 23,
      pt.bg = "white",
      pt.cex = .95,
      pt.lwd = 1.05,
      col = "black",
      horiz = TRUE,
      bty = "n",
      cex = .91,
      lwd = 1.30,
      seg.len = 2.0,
      x.intersp = .65
    )

    mtext(
      "Type I error rate",
      side = 2,
      outer = TRUE,
      line = 1.0,
      at = .60,
      las = 0,
      cex = 1.18
    )

    invisible(NULL)
  }

  plot_figure()

  writeLines(
    c(
      "Null-effect rejection rates in Simulation Study 1 (tau = 0).",
      "Panels correspond to N = 500, 1,000, and 2,000.",
      "Logistic DIF uses Y_T ~ M + G and a two-sided Wald test.",
      "Logistic rejection rates use valid fits as their denominators.",
      "DID results are pooled from raw replications across matching-variable conditions.",
      "The dotted gray horizontal line denotes the nominal significance level."
    ),
    file.path(directory, "figure4_caption.txt")
  )

  cat("\nStudy 1 pooled DID results:\n")
  print(did_pooled, row.names = FALSE, digits = 4)

  cat("\nStudy 1 logistic fit diagnostics:\n")
  print(
    logistic_summary[
      ,
      c(
        keys, "n_valid", "n_failed",
        "n_fits_with_warnings", "n_warnings"
      )
    ],
    row.names = FALSE
  )

  cat("\nStudy 1 outputs complete.\n")
  invisible(NULL)
}


# ==============================================================================
# D. STUDY 2 -- SIMULATION
# ==============================================================================

simulate_study2 <- function(
    replications,
    base_seed,
    directory) {

  prepare_directory(directory)

  if (file.exists(file.path(directory, "simulation.rds"))) {
    message("Study 2 simulation is enabled: existing saved results will be overwritten.")
  }

  mu <- -.5
  mu_assumed <- mu

  forms <- "Smooth"
  eta_values <- c(0, .05, .10, .15)

  calibration <- expand.grid(
    eta_true = eta_values,
    violation = forms,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  calibration <- calibration[
    order(
      match(calibration$violation, forms),
      calibration$eta_true
    ),
    ,
    drop = FALSE
  ]

  rownames(calibration) <- NULL
  calibration$calibration_id <- seq_len(nrow(calibration))
  calibration$eta_assumed <- calibration$eta_true
  calibration$mu <- mu
  calibration$mu_assumed <- mu_assumed

  calibration$parameter <- vapply(
    seq_len(nrow(calibration)),
    function(i) {
      calibrate_parameter(
        eta = calibration$eta_true[i],
        violation = calibration$violation[i]
      )
    },
    numeric(1)
  )

  calibration$parameter_name <- "eta"

  calibration$eta_achieved <- vapply(
    seq_len(nrow(calibration)),
    function(i) {
      calibration$parameter[i]
    },
    numeric(1)
  )

  calibration$id_error <- vapply(
    seq_len(nrow(calibration)),
    function(i) {
      identification_error(
        violation = calibration$violation[i],
        parameter = calibration$parameter[i],
        mu = mu
      )
    },
    numeric(1)
  )

  calibration$b_normal <- normal_bound(
    calibration$eta_assumed,
    calibration$mu_assumed
  )

  calibration$b_df <- distribution_free_bound(
    calibration$eta_assumed
  )
  calibration$absolute_id_error <- abs(calibration$id_error)

  calibration$ratio_normal <- ifelse(
    calibration$b_normal == 0,
    0,
    abs(calibration$id_error) / calibration$b_normal
  )

  calibration$ratio_df <- ifelse(
    calibration$b_df == 0,
    0,
    abs(calibration$id_error) / calibration$b_df
  )

  stopifnot(
    nrow(calibration) == 4L,
    max(abs(calibration$eta_achieved - calibration$eta_true)) < 1e-8,
    all(abs(calibration$id_error) <= calibration$b_normal + 1e-9),
    all(calibration$b_normal <= calibration$b_df + 1e-12)
  )

  design_index <- expand.grid(
    N = c(500L, 1000L, 2000L),
    tau = c(0, -.05, -.10),
    calibration_id = calibration$calibration_id,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  calibration_match <- match(
    design_index$calibration_id,
    calibration$calibration_id
  )

  design <- cbind(
    data.frame(
      condition_id = seq_len(nrow(design_index)),
      N = design_index$N,
      tau = design_index$tau,
      stringsAsFactors = FALSE
    ),
    calibration[calibration_match, , drop = FALSE]
  )

  rownames(design) <- NULL

  gamma_grid <- expand.grid(
    tau = c(0, -.05, -.10),
    calibration_id = calibration$calibration_id,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  gamma_grid$gamma <- vapply(
    seq_len(nrow(gamma_grid)),
    function(i) {
      j <- match(
        gamma_grid$calibration_id[i],
        calibration$calibration_id
      )

      calibrate_gamma(
        tau = gamma_grid$tau[i],
        violation = calibration$violation[j],
        parameter = calibration$parameter[j],
        mu = mu
      )
    },
    numeric(1)
  )

  gamma_grid$tau_achieved <- vapply(
    seq_len(nrow(gamma_grid)),
    function(i) {
      j <- match(
        gamma_grid$calibration_id[i],
        calibration$calibration_id
      )

      tau_from_gamma(
        gamma = gamma_grid$gamma[i],
        violation = calibration$violation[j],
        parameter = calibration$parameter[j],
        mu = mu
      )
    },
    numeric(1)
  )

  gamma_grid$population_DID <- vapply(
    seq_len(nrow(gamma_grid)),
    function(i) {
      j <- match(
        gamma_grid$calibration_id[i],
        calibration$calibration_id
      )

      population_did(
        gamma = gamma_grid$gamma[i],
        violation = calibration$violation[j],
        parameter = calibration$parameter[j],
        mu = mu
      )
    },
    numeric(1)
  )

  gamma_key <- paste(
    gamma_grid$calibration_id,
    gamma_grid$tau,
    sep = "_"
  )

  design_key <- paste(
    design$calibration_id,
    design$tau,
    sep = "_"
  )

  gamma_match <- match(design_key, gamma_key)

  design$gamma <- gamma_grid$gamma[gamma_match]
  design$tau_achieved <- gamma_grid$tau_achieved[gamma_match]
  design$population_DID <- gamma_grid$population_DID[gamma_match]
  design$condition_seed <- as.double(base_seed) +
    1009 * design$condition_id

  stopifnot(
    nrow(design) == 36L,
    !anyNA(gamma_match),
    max(abs(design$tau_achieved - design$tau)) < 1e-8,
    max(abs(
      design$population_DID - design$tau - design$id_error
    )) < 1e-8
  )

  raw_list <- vector("list", nrow(design))

  for (i in seq_len(nrow(design))) {
    d <- design[i, ]
    set_condition_seed(base_seed, d$condition_id)

    did <- numeric(replications)
    did_se <- numeric(replications)
    n0 <- integer(replications)
    n1 <- integer(replications)

    for (r in seq_len(replications)) {
      group <- rbinom(d$N, size = 1L, prob = .5)
      theta <- rnorm(d$N, mean = mu * group, sd = 1)

      y_anchor <- rbinom(
        d$N,
        size = 1L,
        prob = comparison_anchor_irf(theta, d$violation, d$parameter)
      )

      test_probability <- plogis(
        reference_lp(
          theta = theta,
          violation = d$violation,
          parameter = d$parameter
        ) + d$gamma * group
      )

      y_test <- rbinom(
        d$N,
        size = 1L,
        prob = test_probability
      )

      fit <- fit_did_hc3(y_test, y_anchor, group)

      did[r] <- fit["estimate"]
      did_se[r] <- fit["se"]
      n0[r] <- fit["n0"]
      n1[r] <- fit["n1"]
    }

    raw_list[[i]] <- data.frame(
      condition_id = rep(d$condition_id, replications),
      replication = seq_len(replications),
      did = did,
      did_se = did_se,
      n0 = n0,
      n1 = n1,
      stringsAsFactors = FALSE
    )

    cat(
      sprintf(
        "Study 2: condition %d/%d completed (N=%d, tau=%.2f, %s, eta=%.2f).\n",
        i, nrow(design), d$N, d$tau, d$violation, d$eta_true
      )
    )

    flush.console()
  }

  raw <- do.call(rbind, raw_list)
  rownames(raw) <- NULL

  metadata <- make_metadata(
    study = "Study 2",
    replications = replications,
    base_seed = base_seed,
    alpha = ALPHA
  )

  metadata$n_conditions <- nrow(design)
  metadata$eta_values <- eta_values
  metadata$model_id <- STUDY2_MODEL_ID
  metadata$test_reference_irf <- "plogis(1.5 * theta)"
  metadata$test_focal_irf <- "plogis(1.5 * theta + gamma)"
  metadata$anchor_irf <- "eta + (1 - 2 * eta) * plogis(1.5 * theta)"
  metadata$anchor_discrimination <- NULL
  metadata$anchor_difficulty <- NULL
  metadata$tau_values <- c(0, -.05, -.10)
  metadata$sample_sizes <- c(500L, 1000L, 2000L)
  metadata$violation_forms <- forms
  metadata$eta_assumed_equals_eta_true <- TRUE
  metadata$mu_assumed_equals_mu_true <- TRUE
  metadata$external_calibration_uncertainty_included <- FALSE
  metadata$methods <- c(
    "Unadjusted",
    "Normal",
    "Distribution_free"
  )

  simulation <- list(
    raw = raw,
    design = design,
    calibration = calibration,
    gamma_calibration = gamma_grid,
    metadata = metadata
  )

  saveRDS(
    simulation,
    file.path(directory, "simulation.rds")
  )

  write_csv(
    design,
    file.path(directory, "study2_design.csv")
  )

  write_csv(
    calibration,
    file.path(directory, "study2_calibration.csv")
  )

  write_csv(
    gamma_grid,
    file.path(directory, "study2_gamma_calibration.csv")
  )

  save_run_information(metadata, directory)

  cat("Study 2 simulation saved.\n")
  invisible(simulation)
}


# ==============================================================================
# E. STUDY 2 -- OUTPUTS ONLY
# ==============================================================================

outputs_study2 <- function(directory) {
  input_file <- file.path(directory, "simulation.rds")

  if (!file.exists(input_file)) {
    stop("No saved Study 2 simulation. Set RUN_STUDY2_SIMULATION <- TRUE first.")
  }

  simulation <- readRDS(input_file)
  if (!identical(simulation$metadata$model_id, STUDY2_MODEL_ID)) {
    stop("Saved Study 2 results use a different or unversioned DGP. ",
         "Rerun RUN_STUDY2_SIMULATION; legacy results cannot be relabeled.")
  }

  keys <- c(
    "condition_id", "N", "tau", "violation",
    "eta_true", "eta_assumed", "parameter",
    "mu", "mu_assumed", "population_DID",
    "id_error", "absolute_id_error", "b_normal", "b_df"
  )

  required_raw <- c("condition_id", "did", "did_se")

  if (!all(required_raw %in% names(simulation$raw))) {
    stop("Saved Study 2 raw data are missing required columns.")
  }

  if (!all(keys %in% names(simulation$design))) {
    stop("Saved Study 2 design is missing required columns.")
  }

  if (anyDuplicated(simulation$design$condition_id)) {
    stop("Study 2 design contains duplicated condition IDs.")
  }

  if (!all(
    simulation$raw$condition_id %in%
    simulation$design$condition_id
  )) {
    stop("Some raw-data condition IDs are absent from the Study 2 design.")
  }

  results <- merge(
    simulation$raw,
    simulation$design,
    by = "condition_id",
    sort = FALSE
  )

  alpha <- simulation$metadata$alpha

  if (length(alpha) != 1L ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1) {
    stop("Saved Study 2 alpha is invalid.")
  }

  if (any(!is.finite(results$did)) ||
      any(!is.finite(results$did_se)) ||
      any(results$did_se < 0)) {
    stop("Saved Study 2 results contain invalid DID estimates or SEs.")
  }

  groups <- split(results, results$condition_id)

  methods <- c(
    "Unadjusted",
    "Normal",
    "Distribution_free"
  )

  cat(
    "Study 2 outputs use saved replications per condition: ",
    simulation$metadata$replications,
    "\n",
    sep = ""
  )

  saved_counts <- vapply(groups, nrow, integer(1))

  if (length(saved_counts) != nrow(simulation$design)) {
    stop("Some Study 2 design conditions have no saved replications.")
  }

  if (length(simulation$metadata$replications) == 1L &&
      is.finite(simulation$metadata$replications)) {
    if (any(saved_counts != simulation$metadata$replications)) {
      stop("Saved Study 2 replication counts do not match metadata.")
    }
  }

  summary <- do.call(rbind, lapply(groups, function(d) {
    critical <- qnorm(1 - alpha / 2)
    lower <- d$did - critical * d$did_se
    upper <- d$did + critical * d$did_se
    coverage_did <- mean(lower <= d$population_DID[1] &
                           upper >= d$population_DID[1])
    cbind(
      d[1, keys],
      replications = nrow(d),
      coverage_tau_unadjusted = mean(lower <= d$tau[1] & upper >= d$tau[1]),
      coverage_population_DID = coverage_did,
      coverage_population_DID_mcse = binomial_mcse(coverage_did, nrow(d))
    )
  }))

  rownames(summary) <- NULL

  summary <- summary[
    order(summary$condition_id),
    ,
    drop = FALSE
  ]

  intervals <- do.call(rbind, lapply(groups, function(d) {
    s <- sensitivity_summary(
      estimates = d$did,
      standard_errors = d$did_se,
      tau = d$tau[1],
      b_normal = d$b_normal[1],
      b_df = d$b_df[1],
      alpha = alpha
    )

    md <- d[rep(1L, nrow(s)), keys, drop = FALSE]
    rownames(md) <- NULL

    cbind(md, s)
  }))

  rownames(intervals) <- NULL

  intervals <- intervals[
    order(
      intervals$condition_id,
      match(intervals$method, methods)
    ),
    ,
    drop = FALSE
  ]

  stopifnot(
    nrow(summary) == 36L,
    nrow(intervals) == 108L
  )

  # Numerical consistency, not a claim of exact finite-sample
  # nominal coverage.
  for (id in unique(intervals$condition_id)) {
    d <- intervals[
      intervals$condition_id == id,
      ,
      drop = FALSE
    ]

    stopifnot(
      nrow(d) == 3L,
      !anyDuplicated(d$method),
      all(methods %in% d$method)
    )

    d <- d[match(methods, d$method), , drop = FALSE]

    stopifnot(
      all(diff(d$coverage_tau) >= -1e-12),
      all(diff(d$mean_width) >= -1e-12),
      all(diff(d$excludes_zero) <= 1e-12),
      max(abs(
        d$mean_width - d$mean_width[1] - 2 * d$bound
      )) < 1e-10
    )
  }

  table4 <- simulation$design[, keys, drop = FALSE]

  for (method in methods) {
    d <- intervals[
      intervals$method == method,
      ,
      drop = FALSE
    ]

    idx <- match(table4$condition_id, d$condition_id)

    stopifnot(!anyNA(idx))

    for (metric in c(
      "coverage_tau",
      "coverage_mcse",
      "mean_width",
      "excludes_zero"
    )) {
      table4[[paste(method, metric, sep = "_")]] <-
        d[[metric]][idx]
    }
  }

  # Main figure: tau = 0.
  # Other tau values and interval widths are in tables.
  #
  # All display settings belong to outputs_study2(), so the nested
  # plot_figure() can access them without relying on global variables.

  forms <- "Smooth"

  method_lty <- c(
    Unadjusted = 1,
    Normal = 2,
    Distribution_free = 4
  )

  method_pch <- c(
    Unadjusted = 21,
    Normal = 22,
    Distribution_free = 24
  )

  method_labels <- c(
    Unadjusted = "Unadjusted",
    Normal = "Normal bound",
    Distribution_free = "Distribution-free bound"
  )

  method_offsets <- c(
    Unadjusted = 0,
    Normal = -.002,
    Distribution_free = .002
  )

  figure_data <- intervals[
    !is.na(intervals$tau) & intervals$tau == 0,
    ,
    drop = FALSE
  ]

  stopifnot(
    nrow(figure_data) == 36L,
    all(as.character(figure_data$method) %in% methods),
    all(as.character(figure_data$violation) %in% forms),
    all(is.finite(figure_data$eta_true)),
    all(is.finite(figure_data$coverage_tau)),
    all(
      figure_data$coverage_tau >= 0 &
        figure_data$coverage_tau <= 1
    )
  )

  figure_data$series_id <- match(
    as.character(figure_data$method),
    methods
  )

  # Display-only horizontal offsets.
  # Neither eta_true nor coverage_tau is changed.
  figure_data$x_offset <- unname(
    method_offsets[as.character(figure_data$method)]
  )

  figure_data$x_plot <-
    figure_data$eta_true + figure_data$x_offset

  stopifnot(
    all(is.finite(figure_data$x_offset)),
    all(is.finite(figure_data$x_plot)),
    max(abs(
      figure_data$x_offset -
        unname(method_offsets[as.character(figure_data$method)])
    )) < 1e-12
  )

  # Validate every plotted series before opening the graphics device.
  for (form in forms) {
    for (current_N in c(500, 1000, 2000)) {
      for (method in methods) {
        d <- figure_data[
          figure_data$N == current_N &
            figure_data$violation == form &
            figure_data$method == method,
          ,
          drop = FALSE
        ]

        d <- d[order(d$eta_true), , drop = FALSE]

        stopifnot(
          nrow(d) == 4L,
          max(abs(d$eta_true - c(0, .05, .10, .15))) < 1e-12
        )
      }
    }
  }

  # Save numerical outputs before generating the figure.
  write_csv(
    summary,
    file.path(directory, "study2_population_coverage_diagnostic.csv")
  )

  write_csv(
    intervals,
    file.path(directory, "study2_interval_summary_long.csv")
  )

  write_csv(
    table4,
    file.path(directory, "study2_intervals_all_tau_wide.csv")
  )

  # Current appendix: S3 bounds; S4--S6 coverage and mean interval width.
  # The deleted sampling-performance tables are not regenerated.
  bound_table <- simulation$calibration[, c(
    "eta_true", "absolute_id_error", "b_normal", "b_df"
  ), drop = FALSE]
  names(bound_table)[1] <- "eta"
  write_csv(bound_table, file.path(directory, "tableS3_absolute_id_error_bounds.csv"))
  tau_values <- c(0, -.05, -.10)
  for (j in seq_along(tau_values)) {
    tab <- table4[table4$tau == tau_values[j], , drop = FALSE]
    tab <- tab[order(tab$eta_true, tab$N), , drop = FALSE]
    selected <- c("eta_true", "N",
                  paste(methods, "coverage_tau", sep = "_"),
                  paste(methods, "mean_width", sep = "_"))
    tab <- tab[, selected, drop = FALSE]
    names(tab)[1] <- "eta"
    write_csv(tab, file.path(directory, sprintf("tableS%d_intervals.csv", j + 3L)))
  }

  write_csv(
    simulation$calibration,
    file.path(directory, "population_bounds.csv")
  )

  write_csv(
    figure_data,
    file.path(directory, "figure5_data.csv")
  )

  plot_figure <- function() {
    # Shared vertical limits across all three panels.
    # Expand downward automatically if future results fall below .60.
    minimum_coverage <- min(figure_data$coverage_tau)
    y_lower <- .60

    if (minimum_coverage < .60) {
      y_lower <- max(
        0,
        floor((minimum_coverage - .02) / .05) * .05
      )
    }

    y_limits <- c(y_lower, 1.025)

    y_ticks <- sort(unique(c(
      seq(y_lower, .90, by = .10),
      .95,
      1.00
    )))

    y_tick_labels <- sub(
      "^0",
      "",
      sprintf("%.2f", y_ticks)
    )

    pdf(
      file.path(directory, "figure5.pdf"),
      width = 8.5,
      height = 4.3,
      family = "serif",
      useDingbats = FALSE
    )

    on.exit(dev.off(), add = TRUE)

    layout(
      matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE),
      heights = c(3.2, .72)
    )

    par(
      oma = c(0, 2.4, .2, .2),
      cex = 1,
      family = "serif"
    )

    for (current_N in c(500, 1000, 2000)) {
      par(
        mar = c(2.5, 2.8, 2.1, .8),
        mgp = c(1.6, .50, 0),
        tcl = -.23,
        las = 1,
        cex = 1,
        cex.axis = 1.05,
        cex.main = 1.12,
        xaxs = "i",
        yaxs = "i"
      )

      plot(
        NA,
        xlim = c(-.015, .165),
        ylim = y_limits,
        xlab = "",
        ylab = "",
        xaxt = "n",
        yaxt = "n",
        bty = "l"
      )

      title(
        main = bquote(
          italic(N) == .(
            format(
              current_N,
              big.mark = ",",
              trim = TRUE
            )
          )
        ),
        line = .65
      )

      axis(
        side = 1,
        at = c(0, .05, .10, .15),
        labels = c("0", ".05", ".10", ".15")
      )

      axis(
        side = 2,
        at = y_ticks,
        labels = y_tick_labels
      )

      abline(
        h = .95,
        lty = 3,
        lwd = .8,
        col = "grey55"
      )

      # Draw every line first, then every symbol.
      for (method in methods) {
        d <- figure_data[
          figure_data$N == current_N &
            figure_data$method == method,
          ,
          drop = FALSE
        ]

        d <- d[order(d$eta_true), , drop = FALSE]

        lines(
          x = d$x_plot,
          y = d$coverage_tau,
          lty = unname(method_lty[method]),
          lwd = 1.30,
          col = "black"
        )
      }

      for (method in methods) {
        d <- figure_data[
          figure_data$N == current_N &
            figure_data$method == method,
          ,
          drop = FALSE
        ]

        d <- d[order(d$eta_true), , drop = FALSE]

        points(
          x = d$x_plot,
          y = d$coverage_tau,
          pch = unname(method_pch[method]),
          cex = .92,
          lwd = 1.05,
          col = "black",
          bg = "white"
        )
      }
    }

    # Shared bottom strip.
    par(
      mar = c(0, 0, 0, 0),
      cex = 1,
      las = 1
    )

    plot.new()

    text(
      x = .5,
      y = .84,
      labels = expression(
        paste("Maximum absolute IRF difference (", eta, ")")
      ),
      cex = 1.15
    )

    legend(
      x = .5,
      y = .37,
      xjust = .5,
      yjust = .5,
      legend = unname(method_labels[methods]),
      pch = unname(method_pch[methods]),
      lty = unname(method_lty[methods]),
      lwd = 1.30,
      pt.bg = "white",
      pt.cex = .95,
      pt.lwd = 1.05,
      col = "black",
      horiz = TRUE,
      bty = "n",
      cex = 1.02,
      seg.len = 2.0,
      x.intersp = .75
    )

    mtext(
      expression(paste("Coverage of ", tau)),
      side = 2,
      outer = TRUE,
      line = 1.0,
      at = .60,
      las = 0,
      cex = 1.18
    )

    invisible(NULL)
  }

  plot_figure()

  writeLines(
    c(
      "Coverage of the item-bias estimand in Simulation Study 2 (tau = 0).",
      "Panels correspond to N = 500, 1,000, and 2,000.",
      "IRF nonequivalence follows a single smooth path indexed by eta.",
      "Line types and symbols distinguish the unadjusted 95% confidence",
      "interval and sensitivity-adjusted confidence intervals using the normal-distribution",
      "bound or the distribution-free bound.",
      "For an unadjusted interval [L, U], adjusted intervals are",
      "[L - b_N, U + b_N] and [L - 2*eta, U + 2*eta], where",
      "b_N = eta * {4*Phi(abs(mu)/2) - 2}.",
      "Sensitivity inputs equal the data-generating values and are treated as fixed.",
      "Intervals are not clipped to [-1, 1].",
      "The dotted horizontal line denotes .95.",
      "All panels use the same vertical limits.",
      "Normal-bound and distribution-free adjustments are displayed at",
      "eta - .002 and eta + .002, respectively, solely to separate overlaps.",
      "The unadjusted interval is displayed at the true eta.",
      "The true eta levels are 0, .05, .10, and .15.",
      "No vertical offsets are applied.",
      "Results for nonzero tau and interval widths are reported in supplementary tables."
    ),
    file.path(directory, "figure5_caption.txt")
  )

  cat("\nAbsolute identification errors and sensitivity bounds (Table S3):\n")

  print(
    simulation$calibration[
      ,
      c(
        "eta_true", "eta_achieved",
        "absolute_id_error", "b_normal", "b_df", "ratio_normal"
      ),
      drop = FALSE
    ],
    row.names = FALSE,
    digits = 6
  )

  cat("\nSampling-inference diagnostic for the manuscript example:\n")

  print(
    summary[summary$tau == 0 & summary$eta_true == .15 & summary$N == 2000,
            c("tau", "eta_true", "N", "coverage_tau_unadjusted",
              "coverage_population_DID"), drop = FALSE],
    row.names = FALSE,
    digits = 6
  )

  cat("\nMain Figure results: tau = 0\n")

  print(
    figure_data[
      ,
      c(
        "N", "eta_true", "method",
        "coverage_tau", "coverage_mcse", "mean_width"
      ),
      drop = FALSE
    ],
    row.names = FALSE,
    digits = 4
  )

  cat("\nStudy 2 outputs complete. Adjusted coverage may be conservative.\n")
  invisible(NULL)
}


# ==============================================================================
# F. SELF-CONTAINED NUMERICAL CHECKS
# ==============================================================================

check_implementation <- function() {
  # Deterministic sample: these checks do not draw from the simulation RNG.
  group <- rep(c(0L, 1L), c(43L, 57L))
  y_anchor <- rep(c(0, 1, 1, 0, 1), 20)
  y_test <- rep(c(1, 0, 1, 1), 25)
  difference <- y_test - y_anchor

  fast <- fit_did_hc3(y_test, y_anchor, group)

  fit <- lm(difference ~ group)
  X <- model.matrix(fit)
  bread <- solve(crossprod(X))
  adjusted_residual <- residuals(fit) / (1 - hatvalues(fit))

  meat <- crossprod(
    X,
    X * as.numeric(adjusted_residual^2)
  )

  v <- bread %*% meat %*% bread

  stopifnot(
    abs(fast["estimate"] - coef(fit)["group"]) < 1e-12,
    abs(fast["se"] - sqrt(v["group", "group"])) < 1e-12
  )

  grid <- seq(-30, 30, length.out = 100001L)

  for (eta in c(0, .05, .10, .15)) {
    parameter <- calibrate_parameter(eta, "Smooth")

    grid_max <- max(abs(
      reference_irf(grid, "Smooth", parameter) -
        comparison_anchor_irf(grid, "Smooth", parameter)
    ))

    B <- identification_error(
      violation = "Smooth",
      parameter = parameter
    )

    stopifnot(
      abs(parameter - eta) < 1e-12,
      max(abs(reference_irf(grid, "Smooth", parameter) - plogis(1.5 * grid))) < 1e-12,
      max(abs(comparison_anchor_irf(grid, "Smooth", parameter) -
                (eta + (1 - 2 * eta) * plogis(1.5 * grid)))) < 1e-12,
      abs(grid_max - eta) < 1e-6,
      B <= 1e-12,
      abs(B + eta * 0.259964188188494) < 1e-8,
      abs(B) <= normal_bound(eta, -.5) + 1e-9,
      abs(B) <= distribution_free_bound(eta) + 1e-9,
      abs(identification_error(
        violation = "Smooth",
        parameter = parameter,
        mu = 0
      )) < 1e-10
    )

    for (tau in c(0, -.05, -.10)) {
      gamma <- calibrate_gamma(
        tau = tau,
        violation = "Smooth",
        parameter = parameter
      )

      stopifnot(
        abs(gamma - calibrate_gamma(tau, "Difficulty", 0)) < 1e-10,
        abs(
          tau_from_gamma(
            gamma = gamma,
            violation = "Smooth",
            parameter = parameter
          ) - tau
        ) < 1e-8,
        abs(
          population_did(
            gamma = gamma,
            violation = "Smooth",
            parameter = parameter
          ) - tau - B
        ) < 1e-8
      )
    }
  }

  stopifnot(normal_bound(.15, 0) == 0)

  for (mu in c(-1, -.5, .5, 1)) {
    l1 <- integrate(
      function(theta) {
        abs(
          dnorm(theta, mean = mu, sd = 1) -
            dnorm(theta, mean = 0, sd = 1)
        )
      },
      lower = -Inf,
      upper = Inf,
      subdivisions = 1000L,
      rel.tol = 1e-8
    )$value

    stopifnot(
      abs(l1 - (4 * pnorm(abs(mu) / 2) - 2)) < 1e-7
    )
  }

  estimates <- seq(-.15, .10, length.out = 100)
  standard_errors <- rep(.04, 100)

  intervals <- sensitivity_summary(
    estimates = estimates,
    standard_errors = standard_errors,
    tau = 0,
    b_normal = normal_bound(.10, -.5),
    b_df = .20
  )

  stopifnot(
    all(diff(intervals$coverage_tau) >= 0),
    all(diff(intervals$excludes_zero) <= 0),
    all(diff(intervals$mean_width) >= 0),
    max(abs(
      intervals$mean_width -
        intervals$mean_width[1] -
        2 * intervals$bound
    )) < 1e-12
  )

  # Replication-level nesting and exact width expansion.
  critical <- qnorm(.975)
  lower <- estimates - critical * standard_errors
  upper <- estimates + critical * standard_errors

  bounds <- c(0, normal_bound(.10, -.5), .20)

  for (bound in bounds) {
    adjusted_lower <- lower - bound
    adjusted_upper <- upper + bound

    stopifnot(
      all(adjusted_lower <= lower),
      all(adjusted_upper >= upper),
      max(abs(
        (adjusted_upper - adjusted_lower) -
          (upper - lower) -
          2 * bound
      )) < 1e-12
    )
  }

  # Zero sensitivity bounds reproduce the unadjusted interval.
  zero_bound_intervals <- sensitivity_summary(
    estimates = estimates,
    standard_errors = standard_errors,
    tau = 0,
    b_normal = 0,
    b_df = 0
  )

  stopifnot(
    max(abs(
      zero_bound_intervals$coverage_tau -
        zero_bound_intervals$coverage_tau[1]
    )) < 1e-12,
    max(abs(
      zero_bound_intervals$mean_width -
        zero_bound_intervals$mean_width[1]
    )) < 1e-12,
    max(abs(
      zero_bound_intervals$excludes_zero -
        zero_bound_intervals$excludes_zero[1]
    )) < 1e-12
  )

  cat("Numerical implementation checks passed.\n")
  invisible(NULL)
}


# ==============================================================================
# G. EMPIRICAL ILLUSTRATION -- 2023 PIAAC KOREAN SAMPLE
# ==============================================================================

check_empirical_inputs <- function(data_file) {
  # Check prerequisites before starting any potentially lengthy simulations.
  packages <- c("dplyr", "sandwich")
  missing_packages <- packages[!vapply(
    packages, requireNamespace, logical(1), quietly = TRUE
  )]

  if (length(missing_packages) > 0L) {
    stop(
      "Install the required empirical-analysis packages: ",
      paste(missing_packages, collapse = ", "),
      "."
    )
  }

  if (!file.exists(data_file)) {
    stop(
      "PIAAC CSV not found: ", data_file,
      ". Set EMPIRICAL_DATA_FILE to the local Korean public-use CSV, ",
      "or choose RUN_MODE <- \"simulations\" to run without PIAAC data."
    )
  }
  invisible(NULL)
}


run_empirical_illustration <- function(data_file, directory, eta, alpha) {
  # All analyses are unweighted. Confidence intervals and p-values use
  # large-sample normal inference. Sensitivity parameters are treated as fixed.
  stopifnot(
    length(eta) == 1L, is.finite(eta), eta >= 0, eta <= 1,
    length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1
  )

  # G1. Read data and define the analytic sample --------------------------------
  # Download and extract the Korean public-use CSV from the OECD database:
  # https://www.oecd.org/en/data/datasets/PIAAC-2nd-Cycle-Database.html
  # CSV codes do not retain the value labels available in SPSS/SAS files.
  piaac <- read.csv(
    file = data_file,
    sep = ";",
    dec = ".",
    na.strings = c(".", ".n", ".v"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  test_item <- "E320004S"
  anchor_item <- "E320003S"
  pv_vars <- paste0("PVLIT", 1:10)
  reference_code <- 2L
  focal_code <- 4L

  required_variables <- c("AGEG10LFS", test_item, anchor_item, pv_vars)
  missing_variables <- setdiff(required_variables, names(piaac))
  if (length(missing_variables) > 0L) {
    stop("Missing CSV columns: ", paste(missing_variables, collapse = ", "))
  }

  # Retain the two age groups and observed responses to both items.
  dat <- dplyr::filter(
    piaac,
    AGEG10LFS %in% c(reference_code, focal_code),
    !is.na(.data[[test_item]]),
    !is.na(.data[[anchor_item]])
  )
  # Scored responses must be binary. Do not silently recode an unexpected
  # missing-value or response code as an incorrect answer.
  for (item in c(test_item, anchor_item)) {
    if (!is.numeric(dat[[item]]) || any(!dat[[item]] %in% c(0, 1))) {
      stop("Expected observed scores 0 or 1 in ", item,
           ". Check the CSV format and missing-value codes.")
    }
  }
  for (pv in pv_vars) {
    if (!is.numeric(dat[[pv]]) || any(is.infinite(dat[[pv]]))) {
      stop("Expected finite numeric plausible values (or NA) in ", pv, ".")
    }
  }

  dat <- dplyr::mutate(
    dat,
    G = ifelse(AGEG10LFS == reference_code, 0, 1),
    Y_T = ifelse(.data[[test_item]] == 1, 1, 0),
    Y_A = ifelse(.data[[anchor_item]] == 1, 1, 0)
  )

  # Reproduce the manuscript sample: 285 reference and 251 focal respondents.
  stopifnot(
    nrow(dat) == 536L,
    sum(dat$G == 0) == 285L,
    sum(dat$G == 1) == 251L
  )
  z_critical <- qnorm(1 - alpha / 2)

  # G2. Conventional uniform logistic-regression DIF ----------------------------
  # Fit logit Pr(Y_T = 1 | M, G) = beta_0 + beta_1 * M + beta_2 * G.
  # Standardize each plausible value using the reference-group mean and SD.
  fit_logistic_pv <- function(pv) {
    d <- dplyr::filter(dat, !is.na(.data[[pv]]))
    reference_mean <- mean(d[[pv]][d$G == 0])
    reference_sd <- sd(d[[pv]][d$G == 0])
    if (sum(d$G == 0) < 2L || sum(d$G == 1) < 2L ||
        !is.finite(reference_sd) || reference_sd <= 0) {
      stop("Insufficient observations or invalid reference-group SD for ", pv, ".")
    }
    d$M <- (d[[pv]] - reference_mean) / reference_sd

    fit <- glm(Y_T ~ M + G, family = binomial(), data = d)
    if (!isTRUE(fit$converged) || isTRUE(fit$boundary) ||
        !is.finite(coef(fit)["G"]) ||
        !is.finite(vcov(fit)["G", "G"]) || vcov(fit)["G", "G"] <= 0) {
      stop("Invalid logistic-regression fit for ", pv, ".")
    }
    data.frame(
      PV = pv,
      beta = unname(coef(fit)["G"]),
      variance = unname(vcov(fit)["G", "G"])
    )
  }
  pv_results <- dplyr::bind_rows(lapply(pv_vars, fit_logistic_pv))

  # Multiple-imputation combining rules:
  # total variance = mean within-PV variance + (1 + 1/m) * between-PV variance.
  n_pv <- nrow(pv_results)
  logistic_estimate <- mean(pv_results$beta)
  within_variance <- mean(pv_results$variance)
  between_variance <- var(pv_results$beta)
  logistic_se <- sqrt(within_variance + (1 + 1 / n_pv) * between_variance)
  logistic_ci <- logistic_estimate + c(-1, 1) * z_critical * logistic_se
  logistic_p <- 2 * pnorm(
    abs(logistic_estimate / logistic_se), lower.tail = FALSE
  )

  # G3. DID estimation and HC3 inference ----------------------------------------
  # The coefficient on G is the difference between the group means of Y_T - Y_A.
  # Under the identifying assumptions, the population DID equals tau.
  dat$D <- dat$Y_T - dat$Y_A
  did_fit <- lm(D ~ G, data = dat)
  did_vcov <- sandwich::vcovHC(did_fit, type = "HC3")
  did_estimate <- unname(coef(did_fit)["G"])
  did_se <- sqrt(did_vcov["G", "G"])
  did_ci <- did_estimate + c(-1, 1) * z_critical * did_se
  did_p <- 2 * pnorm(abs(did_estimate / did_se), lower.tail = FALSE)

  table1 <- data.frame(
    Method = c("Logistic-regression DIF", "DID"),
    Estimate = c(logistic_estimate, did_estimate),
    SE = c(logistic_se, did_se),
    CI_Lower = c(logistic_ci[1], did_ci[1]),
    CI_Upper = c(logistic_ci[2], did_ci[2]),
    p = c(logistic_p, did_p)
  )

  # G4. Distribution-free sensitivity analysis ---------------------------------
  # Manuscript calibration inputs (OECD, 2013):
  #   E320003: slope = 1.446, difficulty = 0.437 (anchor).
  #   E320004: slope = 1.338, difficulty = 0.399 (test).
  # Under P_j(theta) = plogis(1.7 * a_j * (theta - b_j)), the maximum
  # absolute calibrated IRF difference is approximately 0.0329.
  # EMPIRICAL_ETA = 0.033 is supplied as a fixed sensitivity value.
  # This widens the DID interval; it does not correct the DID point estimate.
  bound_df <- 2 * eta
  ci_df <- c(did_ci[1] - bound_df, did_ci[2] + bound_df)

  # G5. Normal-distribution sensitivity analysis -------------------------------
  # For each PV: (focal mean - reference mean) / reference SD.
  # Average the signed differences first, then take the absolute value.
  # The bound assumes equal-variance normal latent-trait distributions.
  pv_differences <- dplyr::bind_rows(lapply(pv_vars, function(pv) {
    reference <- dat[[pv]][dat$G == 0]
    focal <- dat[[pv]][dat$G == 1]
    mean_difference <- mean(focal, na.rm = TRUE) -
      mean(reference, na.rm = TRUE)

    data.frame(
      PV = pv,
      difference = mean_difference,
      mu = mean_difference / sd(reference, na.rm = TRUE)
    )
  }))

  mu_hat <- mean(pv_differences$mu)
  bound_normal <- eta * (4 * pnorm(abs(mu_hat) / 2) - 2)
  ci_normal <- c(did_ci[1] - bound_normal, did_ci[2] + bound_normal)

  # Retain full precision in computations and saved output; round in the paper.
  sensitivity <- data.frame(
    Specification = c("Distribution-free", "Equal-variance normal"),
    eta = eta,
    mu = c(NA, abs(mu_hat)),
    Bound = c(bound_df, bound_normal),
    CI_Lower = c(ci_df[1], ci_normal[1]),
    CI_Upper = c(ci_df[2], ci_normal[2])
  )

  # G6. Display and save the empirical results ---------------------------------
  cat("\nEmpirical illustration: logistic-regression DIF and DID\n")
  print(table1, digits = 6, row.names = FALSE)
  cat("\nEmpirical illustration: sensitivity-adjusted confidence intervals\n")
  print(sensitivity, digits = 6, row.names = FALSE)

  prepare_directory(directory)
  write_csv(table1, file.path(directory, "table1_empirical_results.csv"))
  write_csv(pv_results, file.path(directory, "illustration_logistic_pv_results.csv"))
  write_csv(pv_differences, file.path(directory, "illustration_pv_literacy_difference.csv"))
  write_csv(sensitivity, file.path(directory, "illustration_sensitivity.csv"))

  metadata <- list(
    study = "2023 PIAAC Korean empirical illustration",
    input_file = normalizePath(data_file, winslash = "/"),
    input_md5 = unname(tools::md5sum(data_file)),
    test_item = test_item,
    anchor_item = anchor_item,
    age_variable = "AGEG10LFS",
    reference_code = reference_code,
    focal_code = focal_code,
    n_reference = sum(dat$G == 0),
    n_focal = sum(dat$G == 1),
    plausible_values = pv_vars,
    weighted = FALSE,
    alpha = alpha,
    eta = eta,
    mu_hat = mu_hat,
    sensitivity_inputs_fixed = TRUE,
    confidence_intervals_clipped = FALSE,
    package_versions = vapply(
      c("dplyr", "sandwich"),
      function(package) as.character(utils::packageVersion(package)),
      character(1)
    ),
    R_version = R.version.string
  )
  save_run_information(metadata, directory)
  print(sessionInfo())

  invisible(list(
    table1 = table1,
    pv_results = pv_results,
    pv_differences = pv_differences,
    sensitivity = sensitivity
  ))
}



# ==============================================================================
# H. EXECUTION -- RUN THE TASKS SELECTED ABOVE
# ==============================================================================

run_selected_tasks <- function() {
  switches <- list(
    RUN_STUDY1_SIMULATION,
    RUN_STUDY1_OUTPUTS,
    RUN_STUDY2_SIMULATION,
    RUN_STUDY2_OUTPUTS,
    RUN_CHECKS,
    RUN_EMPIRICAL
  )

  stopifnot(
    all(vapply(
      switches,
      function(x) {
        is.logical(x) &&
          length(x) == 1L &&
          !is.na(x)
      },
      logical(1)
    ))
  )

  stopifnot(
    is.character(OUTPUT_ROOT),
    length(OUTPUT_ROOT) == 1L,
    !is.na(OUTPUT_ROOT),
    nzchar(OUTPUT_ROOT)
  )

  if (RUN_STUDY1_SIMULATION || RUN_STUDY2_SIMULATION) {
    stopifnot(
      length(N_REP) == 1L,
      is.finite(N_REP),
      N_REP >= 2,
      N_REP == floor(N_REP),
      N_REP <= .Machine$integer.max,
      length(ALPHA) == 1L,
      is.finite(ALPHA),
      ALPHA > 0,
      ALPHA < 1
    )

    for (seed in c(STUDY1_SEED, STUDY2_SEED)) {
      stopifnot(
        is.finite(seed),
        seed >= 1,
        seed == floor(seed),
        seed <= .Machine$integer.max - 1009 * 81
      )
    }

    if (N_REP < 5000L) {
      message("Smoke-test/small-run setting: N_REP = ", N_REP)
    }
  }

  cat("Run mode: ", RUN_MODE, "\n", sep = "")
  cat("Working directory: ", getwd(), "\n", sep = "")

  d1 <- file.path(OUTPUT_ROOT, "study1")
  d2 <- file.path(OUTPUT_ROOT, "study2")

  # Fail early if an output-only request has no saved input.
  if (RUN_STUDY1_OUTPUTS &&
      !RUN_STUDY1_SIMULATION &&
      !file.exists(file.path(d1, "simulation.rds"))) {
    stop(
      "No saved Study 1 results. ",
      "Enable RUN_STUDY1_SIMULATION or check OUTPUT_ROOT."
    )
  }

  if (RUN_STUDY2_OUTPUTS &&
      !RUN_STUDY2_SIMULATION &&
      !file.exists(file.path(d2, "simulation.rds"))) {
    stop(
      "No saved Study 2 results. ",
      "Enable RUN_STUDY2_SIMULATION or check OUTPUT_ROOT."
    )
  }

  if (RUN_EMPIRICAL) {
    check_empirical_inputs(EMPIRICAL_DATA_FILE)
  }

  if (RUN_CHECKS) {
    check_implementation()
  }

  if (RUN_STUDY1_SIMULATION) {
    simulate_study1(
      replications = N_REP,
      base_seed = STUDY1_SEED,
      directory = d1
    )
  }

  if (RUN_STUDY1_OUTPUTS) {
    outputs_study1(d1)
  }

  if (RUN_STUDY2_SIMULATION) {
    simulate_study2(
      replications = N_REP,
      base_seed = STUDY2_SEED,
      directory = d2
    )
  }

  if (RUN_STUDY2_OUTPUTS) {
    outputs_study2(d2)
  }

  if (RUN_EMPIRICAL) {
    run_empirical_illustration(
      data_file = EMPIRICAL_DATA_FILE,
      directory = file.path(OUTPUT_ROOT, "empirical"),
      eta = EMPIRICAL_ETA,
      alpha = ALPHA
    )
  }

  cat("\nAll selected tasks completed.\n")
  invisible(NULL)
}


run_selected_tasks()
