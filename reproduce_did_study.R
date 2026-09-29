# Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach
# Run from the repository folder. See README.md for inputs and output files.
# Sections: settings; shared calculations; simulations; figures/tables; PIAAC; run.

# 1. Settings -----------------------------------------------------------------

RUN_MODE <- "simulations"  # simulations, outputs, empirical, all, checks, smoke
N_REP <- 5000L
STUDY1_SEED <- 2026L
STUDY2_SEED <- 2027L
ALPHA <- .05
OUTPUT_ROOT <- "results"
EMPIRICAL_DATA_FILE <- "prgkorp2.csv"
EMPIRICAL_ETA <- .033
FIGURE_FONT <- "Arial"

# A terminal argument overrides RUN_MODE; Source in RStudio uses the setting above.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Supply one run mode; see README.md.")
  if (length(args) == 1L) RUN_MODE <- args[1L]
}
if (getRversion() < "3.6.0") stop("R >= 3.6.0 is required.")

# 2. Shared calculations --------------------------------------------------

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
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  seed <- as.double(base_seed) + 1009 * as.double(condition_id)

  if (!is.finite(seed) || seed < 1 || 
      seed > .Machine$integer.max) {
    stop("Invalid condition-specific seed.")
  }

  set.seed(as.integer(seed))
  invisible(seed)
}

save_run_information <- function(metadata, directory) {
  writeLines(capture.output(dput(metadata)), file.path(directory, "metadata.txt"))

  writeLines(capture.output(sessionInfo()), file.path(directory, "sessionInfo.txt"))

  invisible(NULL)
}

make_metadata <- function(study, replications, base_seed, alpha) {
  list(study = study, replications = as.integer(replications),
    base_seed = as.integer(base_seed),
    condition_seed_rule = "base_seed + 1009 * condition_id", alpha = alpha,
    rng_kind = "Mersenne-Twister", normal_kind = "Inversion", sample_kind = "Rejection",
    group_probability = .5, mu = -.5, latent_sd = 1, anchor_discrimination = 1.5,
    anchor_difficulty = 0, responses_conditionally_independent = TRUE,
    sensitivity_inputs_fixed = TRUE, confidence_intervals_clipped = FALSE,
    R_version = R.version.string, timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
}

integrate_real <- function(fun) {
  integrate(fun, lower = -Inf, upper = Inf, subdivisions = 1000L, rel.tol = 1e-10,
    abs.tol = 1e-12, stop.on.error = TRUE)$value
}

anchor_irf <- function(theta) {
  plogis(1.5 * theta)
}

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
  if (identical(violation, "Difficulty") && identical(parameter, 0)) {
    return(1.5 * (theta - parameter))
  }
  if (identical(violation, "Smooth") && length(parameter) == 1L && 
      is.finite(parameter) && parameter >= 0 && parameter < .5) {
    return(1.5 * theta)
  }
  stop("Unsupported IRF specification for the manuscript simulations.")
}

reference_irf <- function(theta, violation, parameter) {
  plogis(reference_lp(theta, violation, parameter))
}

calibrate_parameter <- function(eta, violation) {
  if (!identical(violation, "Smooth") || length(eta) != 1L || 
      !is.finite(eta) || eta < 0 || eta >= .5) {
    stop("Study 2 requires the smooth anchor model with 0 <= eta < .5.")
  }
  eta
}

tau_from_gamma <- function(gamma, violation, parameter,
    mu = -.5) {

  integrate_real(function(theta) {
    lp <- reference_lp(theta, violation, parameter)

    (plogis(lp + gamma) - plogis(lp)) *
      dnorm(theta, mean = mu, sd = 1)
  })
}

calibrate_gamma <- function(tau, violation, parameter,
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

  uniroot(objective, interval = c(lower, upper), tol = 1e-12)$root
}

identification_error <- function(violation, parameter,
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

population_did <- function(gamma, violation, parameter,
    mu = -.5) {

  focal_difference <- integrate_real(function(theta) {
    (plogis(reference_lp(theta, violation, parameter) + gamma) -
        comparison_anchor_irf(theta, violation, parameter)) * dnorm(theta, mean = mu, sd = 1)
  })

  reference_difference <- integrate_real(function(theta) {
    (reference_irf(theta, violation, parameter) -
        comparison_anchor_irf(theta, violation, parameter)) * dnorm(theta, mean = 0, sd = 1)
  })

  focal_difference - reference_difference
}

normal_bound <- function(eta, mu) {
  if (any(!is.finite(eta)) || any(eta < 0) || 
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

  if (anyNA(y_test) || anyNA(y_anchor) || anyNA(group) || any(!is.finite(y_test)) || 
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

  c(estimate = estimate, se = sqrt(variance), n0 = n0, n1 = n1)
}

fit_logistic_dif <- function(y_test, matching, group) {
  warning_messages <- character(0)
  error_message <- ""

  fit <- tryCatch(withCallingHandlers(glm(y_test ~ matching + group,
        family = binomial(link = "logit"), control = glm.control(maxit = 50L)),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
)

  invalid_result <- function(reason) {
    list(estimate = NA_real_, se = NA_real_, valid = FALSE, failed = TRUE,
      warning_count = length(warning_messages), warning = length(warning_messages) > 0L,
      reason = reason)
  }

  if (is.null(fit)) {
    return(invalid_result(if (nzchar(error_message)) error_message else "glm_error"))
  }

  if (!isTRUE(fit$converged)) {
    return(invalid_result("nonconvergence"))
  }

  if (isTRUE(fit$boundary)) {
    return(invalid_result("boundary_fit"))
  }

  coefficient_table <- tryCatch(coef(summary(fit)), error = function(e) NULL)

  if (is.null(coefficient_table) || 
      !("group" %in% rownames(coefficient_table))) {
    return(invalid_result("missing_group_coefficient"))
  }

  estimate <- unname(coefficient_table["group", "Estimate"])
  se <- unname(coefficient_table["group", "Std. Error"])

  if (!is.finite(estimate) || !is.finite(se) || 
      se <= 0) {
    return(invalid_result("invalid_estimate_or_se"))
  }

  list(estimate = estimate, se = se, valid = TRUE, failed = FALSE,
    warning_count = length(warning_messages), warning = length(warning_messages) > 0L,
    reason = "")
}

binomial_mcse <- function(probability, n) {
  if (!is.finite(n) || n <= 0) {
    return(rep(NA_real_, length(probability)))
  }

  sqrt(probability * (1 - probability) / n)
}

performance <- function(estimates, standard_errors, truth,
    alpha = .05) {

  if (length(estimates) != length(standard_errors)) {
    stop("Estimate and SE vectors must have equal lengths.")
  }

  if (length(truth) != 1L || !is.finite(truth) || length(alpha) != 1L || !is.finite(alpha) || 
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

  data.frame(n_total = n_total, n_valid = n_valid, n_invalid = n_total - n_valid,
    mean_estimate = mean(estimates), bias = mean(estimates) - truth,
    bias_mcse = empirical_sd / sqrt(n_valid), empirical_sd = empirical_sd, mean_se = mean_se,
    se_sd_ratio = if (empirical_sd > 0) {
      mean_se / empirical_sd
    } else {
      NA_real_
    },
    rmse = sqrt(mean((estimates - truth)^2)), coverage = coverage,
    coverage_mcse = binomial_mcse(coverage, n_valid), mean_width = mean(upper - lower),
    rejection = rejection, rejection_mcse = binomial_mcse(rejection, n_valid),
    excludes_zero_negative = negative, excludes_zero_positive = positive,
    stringsAsFactors = FALSE)
}

sensitivity_summary <- function(estimates, standard_errors, tau, b_normal, b_df,
    alpha = .05) {

  if (length(estimates) != length(standard_errors)) {
    stop("Estimate and SE vectors must have equal lengths.")
  }

  if (length(tau) != 1L || !is.finite(tau) || length(b_normal) != 1L || length(b_df) != 1L || 
      !is.finite(b_normal) || !is.finite(b_df) || b_normal < 0 || b_df < b_normal - 1e-12 || 
      !is.finite(alpha) || alpha <= 0 || 
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

  methods <- c("Unadjusted", "Normal", "Distribution_free")

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

    data.frame(method = methods[j], bound = bounds[j], n_total = n_total, n_valid = n_valid,
      n_invalid = n_total - n_valid, coverage_tau = coverage,
      coverage_mcse = binomial_mcse(coverage, n_valid), mean_width = mean(widths),
      mean_width_mcse = sd(widths) / sqrt(n_valid), excludes_zero = exclusion,
      excludes_zero_mcse = binomial_mcse(exclusion, n_valid),
      excludes_zero_negative = negative,
      excludes_zero_negative_mcse = binomial_mcse(negative, n_valid),
      excludes_zero_positive = positive,
      excludes_zero_positive_mcse = binomial_mcse(positive, n_valid), stringsAsFactors = FALSE
)
  })

  answer <- do.call(rbind, answer)
  rownames(answer) <- NULL
  answer
}

# 3. Study 1: simulation --------------------------------------------------

simulate_study1 <- function(replications, base_seed,
    directory) {

  prepare_directory(directory)

  if (file.exists(file.path(directory, "simulation.rds"))) {
    message("Study 1 simulation is enabled: existing saved results will be overwritten.")
  }

  mu <- -.5

  design <- expand.grid(N = c(500L, 1000L, 2000L), tau = c(0, -.05, -.10),
    delta = c(0, .25, .50), rho_M = c(1, .8, .6), KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)

  design$condition_id <- seq_len(nrow(design))
  design$mu <- mu

  calibration <- data.frame(tau = c(0, -.05, -.10), stringsAsFactors = FALSE)

  calibration$gamma <- vapply(calibration$tau,
    function(tau) {
      calibrate_gamma(tau = tau, violation = "Difficulty", parameter = 0, mu = mu)
    },
    numeric(1))

  calibration$tau_achieved <- vapply(seq_len(nrow(calibration)),
    function(i) {
      tau_from_gamma(gamma = calibration$gamma[i], violation = "Difficulty", parameter = 0,
        mu = mu)
    },
    numeric(1))

  design$gamma <- calibration$gamma[ match(design$tau, calibration$tau) ]

  design$condition_seed <- as.double(base_seed) +
    1009 * design$condition_id

  design <- design[ , c("condition_id", "N", "tau", "delta", "rho_M",
      "mu", "gamma", "condition_seed") ]

  stopifnot(nrow(design) == 81L, max(abs(calibration$tau_achieved - calibration$tau)) < 1e-8
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

      y_anchor <- rbinom(d$N, size = 1L, prob = anchor_irf(theta))

      y_test <- rbinom(d$N, size = 1L, prob = plogis(1.5 * theta + d$gamma * group))

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

    raw_list[[i]] <- data.frame(condition_id = rep(d$condition_id, replications),
      replication = seq_len(replications), did = did, did_se = did_se, n0 = n0, n1 = n1,
      logistic = logistic, logistic_se = logistic_se, logistic_valid = logistic_valid,
      logistic_failed = logistic_failed, logistic_warning = logistic_warning,
      logistic_warning_count = logistic_warning_count,
      logistic_failure_reason = logistic_failure_reason, stringsAsFactors = FALSE)

    cat(sprintf(
        "Study 1: condition %d/%d completed (N=%d, tau=%.2f, delta=%.2f, rho_M=%.2f).\n",
        i, nrow(design), d$N, d$tau, d$delta, d$rho_M))

    flush.console()
  }

  raw <- do.call(rbind, raw_list)
  rownames(raw) <- NULL

  metadata <- make_metadata(study = "Study 1", replications = replications,
    base_seed = base_seed, alpha = ALPHA)

  metadata$n_conditions <- nrow(design)
  metadata$matching_error_variance <- "(1 - rho_M) / rho_M"
  metadata$logistic_model <- "Y_T ~ M + G"
  metadata$logistic_inference <- "Model-based Wald confidence interval"

  simulation <- list(raw = raw, design = design, calibration = calibration,
    metadata = metadata)

  saveRDS(simulation, file.path(directory, "simulation.rds"))

  save_run_information(metadata, directory)

  cat("Study 1 simulation saved.\n")
  invisible(simulation)
}

# 4. Study 2: simulation --------------------------------------------------

simulate_study2 <- function(replications, base_seed,
    directory) {

  prepare_directory(directory)

  if (file.exists(file.path(directory, "simulation.rds"))) {
    message("Study 2 simulation is enabled: existing saved results will be overwritten.")
  }

  mu <- -.5
  mu_assumed <- mu

  forms <- "Smooth"
  eta_values <- c(0, .05, .10, .15)

  calibration <- expand.grid(eta_true = eta_values, violation = forms,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  calibration <- calibration[ order(match(calibration$violation, forms), calibration$eta_true
), , drop = FALSE ]

  rownames(calibration) <- NULL
  calibration$calibration_id <- seq_len(nrow(calibration))
  calibration$eta_assumed <- calibration$eta_true
  calibration$mu <- mu
  calibration$mu_assumed <- mu_assumed

  calibration$parameter <- vapply(seq_len(nrow(calibration)),
    function(i) {
      calibrate_parameter(eta = calibration$eta_true[i], violation = calibration$violation[i]
)
    },
    numeric(1))

  calibration$parameter_name <- "eta"

  calibration$eta_achieved <- vapply(seq_len(nrow(calibration)),
    function(i) {
      calibration$parameter[i]
    },
    numeric(1))

  calibration$id_error <- vapply(seq_len(nrow(calibration)),
    function(i) {
      identification_error(violation = calibration$violation[i],
        parameter = calibration$parameter[i], mu = mu)
    },
    numeric(1))

  calibration$b_normal <- normal_bound(calibration$eta_assumed, calibration$mu_assumed)

  calibration$b_df <- distribution_free_bound(calibration$eta_assumed)
  calibration$absolute_id_error <- abs(calibration$id_error)

  calibration$ratio_normal <- ifelse(calibration$b_normal == 0, 0,
    abs(calibration$id_error) / calibration$b_normal)

  calibration$ratio_df <- ifelse(calibration$b_df == 0, 0,
    abs(calibration$id_error) / calibration$b_df)

  stopifnot(nrow(calibration) == 4L,
    max(abs(calibration$eta_achieved - calibration$eta_true)) < 1e-8,
    all(abs(calibration$id_error) <= calibration$b_normal + 1e-9),
    all(calibration$b_normal <= calibration$b_df + 1e-12))

  design_index <- expand.grid(N = c(500L, 1000L, 2000L), tau = c(0, -.05, -.10),
    calibration_id = calibration$calibration_id, KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)

  calibration_match <- match(design_index$calibration_id, calibration$calibration_id)

  design <- cbind(data.frame(condition_id = seq_len(nrow(design_index)), N = design_index$N,
      tau = design_index$tau, stringsAsFactors = FALSE),
    calibration[calibration_match, , drop = FALSE])

  rownames(design) <- NULL

  gamma_grid <- expand.grid(tau = c(0, -.05, -.10),
    calibration_id = calibration$calibration_id, KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)

  gamma_grid$gamma <- vapply(seq_len(nrow(gamma_grid)),
    function(i) {
      j <- match(gamma_grid$calibration_id[i], calibration$calibration_id)

      calibrate_gamma(tau = gamma_grid$tau[i], violation = calibration$violation[j],
        parameter = calibration$parameter[j], mu = mu)
    },
    numeric(1))

  gamma_grid$tau_achieved <- vapply(seq_len(nrow(gamma_grid)),
    function(i) {
      j <- match(gamma_grid$calibration_id[i], calibration$calibration_id)

      tau_from_gamma(gamma = gamma_grid$gamma[i], violation = calibration$violation[j],
        parameter = calibration$parameter[j], mu = mu)
    },
    numeric(1))

  gamma_grid$population_DID <- vapply(seq_len(nrow(gamma_grid)),
    function(i) {
      j <- match(gamma_grid$calibration_id[i], calibration$calibration_id)

      population_did(gamma = gamma_grid$gamma[i], violation = calibration$violation[j],
        parameter = calibration$parameter[j], mu = mu)
    },
    numeric(1))

  gamma_key <- paste(gamma_grid$calibration_id, gamma_grid$tau, sep = "_")

  design_key <- paste(design$calibration_id, design$tau, sep = "_")

  gamma_match <- match(design_key, gamma_key)

  design$gamma <- gamma_grid$gamma[gamma_match]
  design$tau_achieved <- gamma_grid$tau_achieved[gamma_match]
  design$population_DID <- gamma_grid$population_DID[gamma_match]
  design$condition_seed <- as.double(base_seed) +
    1009 * design$condition_id

  stopifnot(nrow(design) == 36L, !anyNA(gamma_match),
    max(abs(design$tau_achieved - design$tau)) < 1e-8, max(abs(
      design$population_DID - design$tau - design$id_error)) < 1e-8)

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

      y_anchor <- rbinom(d$N, size = 1L,
        prob = comparison_anchor_irf(theta, d$violation, d$parameter))

      test_probability <- plogis(reference_lp(theta = theta, violation = d$violation,
          parameter = d$parameter) + d$gamma * group)

      y_test <- rbinom(d$N, size = 1L, prob = test_probability)

      fit <- fit_did_hc3(y_test, y_anchor, group)

      did[r] <- fit["estimate"]
      did_se[r] <- fit["se"]
      n0[r] <- fit["n0"]
      n1[r] <- fit["n1"]
    }

    raw_list[[i]] <- data.frame(condition_id = rep(d$condition_id, replications),
      replication = seq_len(replications), did = did, did_se = did_se, n0 = n0, n1 = n1,
      stringsAsFactors = FALSE)

    cat(sprintf("Study 2: condition %d/%d completed (N=%d, tau=%.2f, %s, eta=%.2f).\n",
        i, nrow(design), d$N, d$tau, d$violation, d$eta_true))

    flush.console()
  }

  raw <- do.call(rbind, raw_list)
  rownames(raw) <- NULL

  metadata <- make_metadata(study = "Study 2", replications = replications,
    base_seed = base_seed, alpha = ALPHA)

  metadata$n_conditions <- nrow(design)
  metadata$eta_values <- eta_values
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
  metadata$methods <- c("Unadjusted", "Normal", "Distribution_free")

  simulation <- list(raw = raw, design = design, calibration = calibration,
    gamma_calibration = gamma_grid, metadata = metadata)

  saveRDS(simulation, file.path(directory, "simulation.rds"))

  save_run_information(metadata, directory)

  cat("Study 2 simulation saved.\n")
  invisible(simulation)
}

# 5. Figures 4 and 5 ------------------------------------------------------

plot_rejection_grid <- function(data, directory, study, alpha = .05) {
  taus <- c(0, -.05, -.10); Ns <- c(500, 1000, 2000)
  near <- function(x, y) abs(x-y) < 1e-8
  rates <- function(d, variable = "rejection") {
    if (nrow(d) != 3L) stop("Expected three effect sizes in each series.")
    ii <- vapply(taus, function(t) {
      z <- which(near(d$tau, t)); if (length(z) != 1L) stop("Missing/duplicated tau."); z
    }, integer(1))
    y <- d[[variable]][ii]
    if (length(y) != 3L || any(!is.finite(y)) || any(y<0 | y>1)) stop("Invalid rates.")
    y
  }
  open_device <- function(path, height) {
    # Native macOS PDF avoids X11/Cairo and improves Greek-font portability.
    if (identical(Sys.info()[["sysname"]], "Darwin") && isTRUE(capabilities("aqua"))) {
      grDevices::quartz(type = "pdf", file = path, width = 9.4, height = height,
                       family = FIGURE_FONT, pointsize = 12)
    } else if (isTRUE(capabilities("cairo"))) {
      grDevices::cairo_pdf(path, width = 9.4, height = height, family = FIGURE_FONT, pointsize = 12)
    } else {
      stop("PDF output requires native macOS Quartz or Cairo support in R.")
    }
  }
  blank <- function() {
    # Fix coordinates explicitly: axes() must not change subsequent strip alignment.
    par(mar = c(0, 0, 0, 0), cex = 1, xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  }
  axes <- function(title = NULL, row_label = NULL) {
    par(mar = c(2.3, 2.8, 2.5, .8), mgp = c(1.7, .5, 0), tcl = -.2, las = 1, xaxs = "i", yaxs = "i", cex = 1)
    plot(NA, xlim = c(.005, -.105), ylim = c(-.025, 1.025), xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "l")
    abline(h = c(.25, .5, .75, 1), col = "grey92", lwd = .6)
    abline(h = alpha, col = "grey60", lty = 3, lwd = .8)
    axis(1, at = taus, labels = c("0", "-0.05", "-0.10"), cex.axis = .95)
    axis(2, at = seq(0, 1, .2), labels = c("0", ".2", ".4", ".6", ".8", "1.0"), cex.axis = .95)
    if (!is.null(title)) mtext(title, side = 3, line = .8, cex = 1.05, las = 1)
    if (!is.null(row_label)) mtext(row_label, side = 2, line = 3.1, las = 1, adj = 1, cex = 1.05, xpd = NA)
  }
  ylabel <- function(at = .55) mtext("Rejection rate", side = 2, outer = TRUE,
                                  line = 2.5, las = 0, at = at, cex = 1.2)
  tau_strip <- function() {blank(); text(.5, .5, expression(tau), cex = 1.15)}
  draw <- function(ys, cols, ltys, pchs, filled, widths, order = seq_along(ys)) {
    for (j in order) lines(taus, ys[[j]], col = cols[j], lty = ltys[j], lwd = widths[j])
    for (j in order) points(taus, ys[[j]], col = cols[j], pch = pchs[j],
                          bg = if (filled[j]) cols[j] else "white", cex = .95, lwd = 1.1)
  }
  verify_pdf <- function(path) {
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) stop("No PDF created: ", path)
    message("Figure saved: ", normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  path <- file.path(directory, paste0("figure", study+3L, ".pdf"))
  plot4 <- function(path) {
    # Validate the full crossed design before drawing.
    for (dd in c(0, .25, .5)) for (r in c(1, .8, .6)) for (nn in Ns) {
      z <- data[near(data$delta, dd) & near(data$rho_M, r) & data$N == nn, , drop = FALSE]
      rates(z); rates(z, "did_rejection")
    }
    open_device(path, 8.85)
    on.exit(dev.off(), add = TRUE)
    par(oma = c(.4, 4.8, .4, .3), family = FIGURE_FONT)
    # For each block: heading, three panels, tau, condition legend.
    layout(rbind(c(1, 1, 1), c(2, 3, 4), c(5, 5, 5), c(6, 6, 6),
                 c(0, 0, 0), # Small empty gap between blocks A and B.
                 c(7, 7, 7), c(8, 9, 10), c(11, 11, 11), c(12, 12, 12)),
           heights = c(.32, 2.3, .24, .48, .10, .32, 2.3, .24, .48))
    cols <- c("black", "grey35", "grey45", "grey55")
    ltys <- c(1, 2, 4, 5); pchs <- c(21, 21, 24, 22); widths <- c(1.9, 1.3, 1.3, 1.3)
    for (block in 1:2) {
      blank()
      # Reset cex in blank(): layout() otherwise shrinks the first heading.
      heading <- if (block == 1) "A" else "B"
      text(0, .5, heading, adj = c(0, .5), cex = 1.2, font = 2)
      for (nn in Ns) {
        axes(bquote(italic(N) == .(format(nn, big.mark = ",", trim = TRUE))))
        oracle <- data[data$N == nn & near(data$delta, 0) & near(data$rho_M, 1), , drop = FALSE]
        ys <- list(rates(oracle, "did_rejection"))
        for (k in 1:3) {
          dd <- if (block == 1) c(0, .25, .5)[k] else 0
          r <- if (block == 1) 1 else c(1, .8, .6)[k]
          z <- data[data$N == nn & near(data$delta, dd) & near(data$rho_M, r), , drop = FALSE]
          ys[[k+1]] <- rates(z)
        }
        draw(ys, cols, ltys, pchs, c(TRUE, FALSE, FALSE, FALSE), widths, order = c(2, 3, 4, 1))
      }
      tau_strip(); blank()
      labs <- if (block == 1) expression(DID, paste("Logistic: ", delta == 0),
                    paste("Logistic: ", delta == .25), paste("Logistic: ", delta == .50)) else
               expression(DID, paste("Logistic: ", rho == 1.00),
                    paste("Logistic: ", rho == .80), paste("Logistic: ", rho == .60))
      legend("center", legend = labs, horiz = TRUE, bty = "n", cex = 1.0,
             col = cols, lty = ltys, lwd = widths, pch = pchs, pt.bg = c("black", rep("white", 3)), seg.len = 1.8)
    }
    ylabel(.55)
  }
  plot5 <- function(path) {
    methods <- c("Unadjusted", "Normal", "Distribution_free");etas <- c(.05, .10, .15)
    for (e in etas) for (nn in Ns) for (m in methods) rates(data[near(data$eta_true, e) & data$N == nn & data$method == m, , drop = FALSE])
    open_device(path, 9.2);on.exit(dev.off(), add = TRUE)
    layout(rbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9), c(10, 10, 10), c(11, 11, 11)),
           heights = c(2.3, 2.3, 2.3, .3, .55))
    par(oma = c(.4, 6.6, .4, .3), family = FIGURE_FONT)
    cols <- c("grey45", "black", "grey30");ltys <- c(2, 1, 4);pchs <- c(21, 21, 24);widths <- c(1.3, 1.9, 1.3)
    for (rr in 1:3) for (cc in 1:3) {
      axes(if (rr == 1) bquote(italic(N) == .(format(Ns[cc], big.mark = ",", trim = TRUE))) else NULL,
           if (cc == 1) bquote(eta == .(sprintf("%.2f", etas[rr]))) else NULL)
      ys <- lapply(methods, function(m) rates(data[near(data$eta_true, etas[rr]) & data$N == Ns[cc] & data$method == m, , drop = FALSE]))
      draw(ys, cols, ltys, pchs, c(FALSE, TRUE, FALSE), widths, c(1, 3, 2))
    }
    tau_strip();blank()
    legend("center", legend = c("Unadjusted", "Normal bound", "Distribution-free bound"),
           horiz = TRUE, bty = "n", cex = 1.05, col = cols, lty = ltys, lwd = widths, pch = pchs, pt.bg = c("white", "black", "white"))
    mtext("Rejection rate", side = 2, outer = TRUE, line = 4.8, las = 0, at = .55, cex = 1.2)
  }
  if (study == 1L) {
    plot4(path)
  } else if (study == 2L) {
    plot5(path)
  } else {
    stop("Unknown study.")
  }
  verify_pdf(path)
  invisible(path)
}

# 6. Appendix tables and figure data --------------------------------------

# CUP table markup used in the manuscript; requires its template and booktabs.
write_appendix_tables <- function(did, logistic, directory) {
  decimal <- function(x, digits = 3L) {
    sub("^(-?)0\\.", "\\1.", sprintf(paste0("%.", digits, "f"), x))
  }
  row <- function(values) paste0(paste(values, collapse = " & "), " \\\\")
  sample_label <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
  start <- function(caption, label, columns) c("\\begin{table}[!htbp]", "\\tabcolsep=0pt",
    paste0("\\TBL{\\caption{", caption, "%"), paste0("\\label{", label, "}}}"),
    "{\\begin{fntable}",
    paste0("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}", columns, "@{}}"),
    "\\toprule")
  finish <- function(note) c("\\botrule", "\\end{tabular*}",
    paste0("\\footnotetext[]{\\textit{Note:} ", note, "}"), "\\end{fntable}}", "\\end{table}"
)

  did <- did[order(did$N, -did$tau), ]
  lines <- start("DID Performance in Study~1", "tab:sim1_did", "rrrrrrr")
  lines <- c(lines, row(c("\\TCH{$N$}", "\\TCH{$\\tau$}", "\\TCH{Bias}", "\\TCH{RMSE}",
    "\\TCH{SE/SD}", "\\TCH{Coverage}", "\\TCH{Rejection}")), "\\midrule")
  for (i in seq_len(nrow(did))) {
    d <- did[i, ]
    new_N <- i == 1L || d$N != did$N[i - 1L]
    if (new_N && i > 1L) lines <- c(lines, "\\addlinespace")
    bias <- round(d$bias, 4)
    if (bias == 0) bias <- 0  # Avoid printing negative zero.
    lines <- c(lines, row(c(if (new_N) sample_label(d$N) else "", decimal(d$tau, 2),
      decimal(bias, 4), decimal(d$rmse), sprintf("%.3f", d$se_sd_ratio),
      decimal(d$coverage, 4), decimal(d$rejection, 4))))
  }
  counts <- unique(did$n_valid)
  count_note <- if (length(counts) == 1L) {
    paste0("Each row pools ", sample_label(counts), " valid replications")
  } else {
    "Each row pools valid replications"
  }
  note <- paste0(
    count_note, " across nine matching-variable conditions, which affect neither ",
    "the DID estimator nor the distribution of its inputs. ",
    "Bias and RMSE are relative to $\\tau$. SE/SD is the ratio of the mean ",
    "estimated standard error to the empirical standard deviation of the pooled estimates. ",
    "Coverage refers to nominal 95\\% confidence intervals for $\\tau$. ",
    "Rejection denotes Type~I error when $\\tau=0$ and power when $\\tau\\neq0$, ",
    "at the two-sided .05 level.")
  writeLines(c(lines, finish(note)), file.path(directory, "table_study1_did.tex"))

  lines <- start("Rejection Rates for Logistic-Regression DIF in Study~1",
                 "tab:sim1_logistic", "rrccccccccc")
  lines <- c(lines, row(c("", "", "\\multicolumn{3}{c}{$\\tau=0$}",
          "\\multicolumn{3}{c}{$\\tau=-.05$}", "\\multicolumn{3}{c}{$\\tau=-.10$}")),
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
    row(c("\\TCH{$N$}", "\\TCH{$\\delta$}", rep("\\multicolumn{3}{c}{$\\rho$}", 3))),
    row(c("", "", rep(c("1.00", ".80", ".60"), 3))), "\\midrule")
  for (n in c(500, 1000, 2000)) {
    if (n != 500) lines <- c(lines, "\\addlinespace")
    for (delta in c(0, .25, .50)) {
      values <- numeric(0)
      for (tau in c(0, -.05, -.10)) {
        for (rho in c(1, .8, .6)) {
          d <- logistic[logistic$N == n & logistic$delta == delta & 
                          logistic$tau == tau & logistic$rho == rho, ]
          if (nrow(d) != 1L) stop("Missing or duplicated logistic table condition.")
          values <- c(values, d$rejection)
        }
      }
      lines <- c(lines, row(c(
        if (delta == 0) sample_label(n) else "", decimal(delta, 2), decimal(values))))
    }
  }
  counts <- unique(logistic$n_valid)
  note <- if (length(counts) == 1L) {
    paste0("Each cell is based on ", sample_label(counts), " valid replications. ")
  } else {
    "Rejection rates use valid fits as denominators; counts are provided in the CSV. "
  }
  note <- paste0(note, "Rejection denotes Type~I error when $\\tau=0$ and power ",
                 "when $\\tau\\neq0$, at the two-sided .05 level.")
  writeLines(c(lines, finish(note)), file.path(directory, "table_study1_logistic.tex"))
}

outputs_study1 <- function(directory) {
  input_file <- file.path(directory, "simulation.rds")

  if (!file.exists(input_file)) {
    stop("No saved Study 1 simulation. Run mode 'simulations' first.")
  }

  simulation <- readRDS(input_file)

  required_raw <- c("condition_id", "did", "did_se", "logistic", "logistic_se",
    "logistic_valid", "logistic_failed", "logistic_warning", "logistic_warning_count")

  required_design <- c("condition_id", "N", "tau", "delta", "rho_M", "gamma")

  if (!all(required_raw %in% names(simulation$raw)) || 
      !all(required_design %in% names(simulation$design))) {
    stop("Saved Study 1 data have a different schema. ",
      "Use saved results from the current simulation design.")
  }

  results <- merge(simulation$raw, simulation$design, by = "condition_id", sort = FALSE)

  alpha <- simulation$metadata$alpha
  groups <- split(results, results$condition_id)

  keys <- c("condition_id", "N", "tau", "delta", "rho_M", "gamma")

  cat("Study 1 outputs use saved replications per condition: ",
    simulation$metadata$replications, "\n", sep = "")

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

    cbind(d[1, keys], data.frame(n_total = n_total, n_valid = n_valid,
        n_failed = n_total - n_valid, n_reported_failures = sum(d$logistic_failed),
        n_fits_with_warnings = sum(d$logistic_warning),
        n_warnings = sum(d$logistic_warning_count), mean_coefficient = mean_coefficient,
        coefficient_sd = coefficient_sd, mean_se = mean_se, rejection = rejection,
        rejection_mcse = binomial_mcse(rejection, n_valid), excludes_zero_negative = negative,
        excludes_zero_positive = positive, stringsAsFactors = FALSE))
  }))

  rownames(logistic_summary) <- NULL

  # Pool raw DID replications, rather than averaging condition-level
  # nonlinear statistics such as empirical SD or RMSE.
  pooled_keys <- unique(results[, c("N", "tau")])
  pooled_keys <- pooled_keys[ order(pooled_keys$N, -pooled_keys$tau), , drop = FALSE ]

  did_pooled <- do.call(rbind,
    lapply(seq_len(nrow(pooled_keys)), function(i) {
      current_N <- pooled_keys$N[i]
      current_tau <- pooled_keys$tau[i]

      d <- results[ results$N == current_N & results$tau == current_tau, , drop = FALSE ]

      cbind(pooled_keys[i, , drop = FALSE], data.frame(
          n_conditions_pooled = length(unique(d$condition_id))), performance(
          estimates = d$did, standard_errors = d$did_se, truth = current_tau, alpha = alpha)
)
    })
)

  rownames(did_pooled) <- NULL

  figure_data <- logistic_summary
  key <- function(N, tau) paste(N, sprintf("%.4f", tau), sep = ":")
  idx <- match(key(figure_data$N, figure_data$tau), key(did_pooled$N, did_pooled$tau))
  stopifnot(!anyNA(idx))
  figure_data$did_rejection <- did_pooled$rejection[idx]
  figure_data$did_rejection_mcse <- did_pooled$rejection_mcse[idx]
  figure_data$metric <- ifelse(figure_data$tau == 0, "Type I error", "Power")

  did_table <- did_pooled[, c(
    "N", "tau", "bias", "rmse", "se_sd_ratio", "coverage", "rejection", "n_valid")]
  logistic_table <- logistic_summary[, c(
    "N", "tau", "delta", "rho_M", "rejection", "n_valid", "n_failed", "n_fits_with_warnings"
)]
  names(logistic_table)[names(logistic_table) == "rho_M"] <- "rho"
  logistic_table <- logistic_table[order(
    logistic_table$N, logistic_table$delta, -logistic_table$tau, -logistic_table$rho), ]
  write_csv(did_table, file.path(directory, "table_study1_did.csv"))
  write_csv(logistic_table, file.path(directory, "table_study1_logistic.csv"))
  write_appendix_tables(did_table, logistic_table, directory)
  write_csv(figure_data, file.path(directory, "figure4_data.csv"))
  plot_rejection_grid(figure_data, directory, study = 1L, alpha = alpha)
  message("Study 1: Figure 4 and Tables S1-S2 saved in ", directory)
  invisible(list(did = did_table, logistic = logistic_table))
}

outputs_study2 <- function(directory) {
  input <- file.path(directory, "simulation.rds")
  if (!file.exists(input)) stop("No saved Study 2 results. Run mode 'simulations' first.")
  simulation <- readRDS(input)
  # Confirm that the saved results use the Study 2 IRFs.
  expected_irfs <- list(
    test_reference_irf = "plogis(1.5 * theta)",
    test_focal_irf = "plogis(1.5 * theta + gamma)",
    anchor_irf = "eta + (1 - 2 * eta) * plogis(1.5 * theta)"
  )
  if (!identical(simulation$metadata[names(expected_irfs)], expected_irfs)) {
    stop("Saved Study 2 results do not match the specified IRFs. Run the simulations first.")
  }
  keys <- c("condition_id", "N", "tau", "eta_true", "b_normal", "b_df")
  if (!all(keys %in% names(simulation$design)) || 
      !all(c("condition_id", "did", "did_se") %in% names(simulation$raw)) || 
      anyDuplicated(simulation$design$condition_id) || 
      !all(simulation$raw$condition_id %in% simulation$design$condition_id)) {
    stop("Invalid saved Study 2 data.")
  }
  results <- merge(simulation$raw, simulation$design, by = "condition_id", sort = FALSE)
  groups <- split(results, results$condition_id)
  if (length(groups) != 36L || 
      any(vapply(groups, nrow, integer(1)) != simulation$metadata$replications) || 
      any(!is.finite(results$did)) || any(!is.finite(results$did_se)) || 
      any(results$did_se < 0)) stop("Incomplete or invalid Study 2 replications.")
  alpha <- simulation$metadata$alpha
  figure_data <- do.call(rbind, lapply(groups, function(d) {
    s <- sensitivity_summary(d$did, d$did_se, d$tau[1], d$b_normal[1], d$b_df[1], alpha)
    data.frame(condition_id = d$condition_id[1], N = d$N[1], tau = d$tau[1],
      eta_true = d$eta_true[1], method = s$method, bound = s$bound,
      n_valid = s$n_valid, rejection = s$excludes_zero, rejection_mcse = s$excludes_zero_mcse,
      metric = if (d$tau[1] == 0) "Type I error" else "Power")
  }))
  rownames(figure_data) <- NULL
  stopifnot(nrow(figure_data) == 108L)
  for (id in unique(figure_data$condition_id)) {
    d <- figure_data[figure_data$condition_id == id, ]
    stopifnot(all(diff(d$rejection) <= 1e-12))
    if (d$eta_true[1] == 0) stopifnot(length(unique(d$rejection)) == 1L)
  }
  write_csv(figure_data, file.path(directory, "figure5_data.csv"))
  plot_rejection_grid(figure_data, directory, study = 2L, alpha = alpha)
  message("Study 2: Figure 5 saved in ", directory)
  invisible(figure_data)
}

# 7. Numerical implementation checks --------------------------------------

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

  meat <- crossprod(X, X * as.numeric(adjusted_residual^2))

  v <- bread %*% meat %*% bread

  stopifnot(abs(fast["estimate"] - coef(fit)["group"]) < 1e-12,
    abs(fast["se"] - sqrt(v["group", "group"])) < 1e-12)

  grid <- seq(-30, 30, length.out = 100001L)

  for (eta in c(0, .05, .10, .15)) {
    parameter <- calibrate_parameter(eta, "Smooth")

    grid_max <- max(abs(reference_irf(grid, "Smooth", parameter) -
        comparison_anchor_irf(grid, "Smooth", parameter)))

    B <- identification_error(violation = "Smooth", parameter = parameter)

    stopifnot(abs(parameter - eta) < 1e-12,
      max(abs(reference_irf(grid, "Smooth", parameter) - plogis(1.5 * grid))) < 1e-12,
      max(abs(comparison_anchor_irf(grid, "Smooth", parameter) -
                (eta + (1 - 2 * eta) * plogis(1.5 * grid)))) < 1e-12,
      abs(grid_max - eta) < 1e-6, B <= 1e-12, abs(B + eta * 0.259964188188494) < 1e-8,
      abs(B) <= normal_bound(eta, -.5) + 1e-9, abs(B) <= distribution_free_bound(eta) + 1e-9,
      abs(identification_error(violation = "Smooth", parameter = parameter, mu = 0)) < 1e-10
)

    for (tau in c(0, -.05, -.10)) {
      gamma <- calibrate_gamma(tau = tau, violation = "Smooth", parameter = parameter)

      stopifnot(abs(gamma - calibrate_gamma(tau, "Difficulty", 0)) < 1e-10, abs(
          tau_from_gamma(gamma = gamma, violation = "Smooth", parameter = parameter) - tau
) < 1e-8, abs(population_did(gamma = gamma, violation = "Smooth",
            parameter = parameter) - tau - B) < 1e-8)
    }
  }

  stopifnot(normal_bound(.15, 0) == 0)

  for (mu in c(-1, -.5, .5, 1)) {
    l1 <- integrate(
      function(theta) {
        abs(dnorm(theta, mean = mu, sd = 1) - dnorm(theta, mean = 0, sd = 1))
      },
      lower = -Inf, upper = Inf, subdivisions = 1000L, rel.tol = 1e-8)$value

    stopifnot(abs(l1 - (4 * pnorm(abs(mu) / 2) - 2)) < 1e-7)
  }

  estimates <- seq(-.15, .10, length.out = 100)
  standard_errors <- rep(.04, 100)

  intervals <- sensitivity_summary(estimates = estimates, standard_errors = standard_errors,
    tau = 0, b_normal = normal_bound(.10, -.5), b_df = .20)

  stopifnot(all(diff(intervals$coverage_tau) >= 0), all(diff(intervals$excludes_zero) <= 0),
    all(diff(intervals$mean_width) >= 0), max(abs(intervals$mean_width -
        intervals$mean_width[1] - 2 * intervals$bound)) < 1e-12)

  # Replication-level nesting and exact width expansion.
  critical <- qnorm(.975)
  lower <- estimates - critical * standard_errors
  upper <- estimates + critical * standard_errors

  bounds <- c(0, normal_bound(.10, -.5), .20)

  for (bound in bounds) {
    adjusted_lower <- lower - bound
    adjusted_upper <- upper + bound

    stopifnot(all(adjusted_lower <= lower), all(adjusted_upper >= upper), max(abs(
        (adjusted_upper - adjusted_lower) - (upper - lower) - 2 * bound)) < 1e-12)
  }

  # Zero sensitivity bounds reproduce the unadjusted interval.
  zero_bound_intervals <- sensitivity_summary(estimates = estimates,
    standard_errors = standard_errors, tau = 0, b_normal = 0, b_df = 0)

  stopifnot(max(abs(zero_bound_intervals$coverage_tau - zero_bound_intervals$coverage_tau[1]
)) < 1e-12, max(abs(zero_bound_intervals$mean_width - zero_bound_intervals$mean_width[1]
)) < 1e-12, max(abs(zero_bound_intervals$excludes_zero -
        zero_bound_intervals$excludes_zero[1])) < 1e-12)

  cat("Numerical implementation checks passed.\n")
  invisible(NULL)
}

# 8. PIAAC empirical illustration (optional) ------------------------------

check_empirical_inputs <- function(data_file) {
  # Check prerequisites before starting any potentially lengthy simulations.
  packages <- c("dplyr", "sandwich")
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE
)]

  if (length(missing_packages) > 0L) {
    stop("Install the required empirical-analysis packages: ",
      paste(missing_packages, collapse = ", "), ".")
  }

  if (!file.exists(data_file)) {
    stop("PIAAC CSV not found: ", data_file,
      ". Set EMPIRICAL_DATA_FILE to the local Korean public-use CSV, ",
      "or choose RUN_MODE <- \"simulations\" to run without PIAAC data.")
  }
  invisible(NULL)
}

run_empirical_illustration <- function(data_file, directory, eta, alpha) {
  # All analyses are unweighted. Confidence intervals and p-values use
  # large-sample normal inference. Sensitivity parameters are treated as fixed.
  stopifnot(length(eta) == 1L, is.finite(eta), eta >= 0, eta <= 1,
    length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1)

  # G1. Read data and define the analytic sample --------------------------------
  # Download and extract the Korean public-use CSV from the OECD database:
  # https://www.oecd.org/en/data/datasets/PIAAC-2nd-Cycle-Database.html
  # CSV codes do not retain the value labels available in SPSS/SAS files.
  piaac <- read.csv(file = data_file, sep = ";", dec = ".", na.strings = c(".", ".n", ".v"),
    check.names = FALSE, stringsAsFactors = FALSE)

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
  dat <- dplyr::filter(piaac, AGEG10LFS %in% c(reference_code, focal_code),
    !is.na(.data[[test_item]]), !is.na(.data[[anchor_item]]))
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

  dat <- dplyr::mutate(dat, G = ifelse(AGEG10LFS == reference_code, 0, 1),
    Y_T = ifelse(.data[[test_item]] == 1, 1, 0), Y_A = ifelse(.data[[anchor_item]] == 1, 1, 0)
)

  # Reproduce the manuscript sample: 285 reference and 251 focal respondents.
  stopifnot(nrow(dat) == 536L, sum(dat$G == 0) == 285L, sum(dat$G == 1) == 251L)
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
    if (!isTRUE(fit$converged) || isTRUE(fit$boundary) || !is.finite(coef(fit)["G"]) || 
        !is.finite(vcov(fit)["G", "G"]) || vcov(fit)["G", "G"] <= 0) {
      stop("Invalid logistic-regression fit for ", pv, ".")
    }
    data.frame(PV = pv, beta = unname(coef(fit)["G"]), variance = unname(vcov(fit)["G", "G"])
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
  logistic_p <- 2 * pnorm(abs(logistic_estimate / logistic_se), lower.tail = FALSE)

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

  table1 <- data.frame(Method = c("Logistic-regression DIF", "DID"),
    Estimate = c(logistic_estimate, did_estimate), SE = c(logistic_se, did_se),
    CI_Lower = c(logistic_ci[1], did_ci[1]), CI_Upper = c(logistic_ci[2], did_ci[2]),
    p = c(logistic_p, did_p))

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

    data.frame(PV = pv, difference = mean_difference,
      mu = mean_difference / sd(reference, na.rm = TRUE))
  }))

  mu_hat <- mean(pv_differences$mu)
  bound_normal <- eta * (4 * pnorm(abs(mu_hat) / 2) - 2)
  ci_normal <- c(did_ci[1] - bound_normal, did_ci[2] + bound_normal)

  # Retain full precision in computations and saved output; round in the paper.
  sensitivity <- data.frame(Specification = c("Distribution-free", "Equal-variance normal"),
    eta = eta, mu = c(NA, abs(mu_hat)), Bound = c(bound_df, bound_normal),
    CI_Lower = c(ci_df[1], ci_normal[1]), CI_Upper = c(ci_df[2], ci_normal[2]))

  # G6. Display and save the empirical results ---------------------------------
  cat("\nEmpirical illustration: logistic-regression DIF and DID\n")
  print(table1, digits = 6, row.names = FALSE)
  cat("\nEmpirical illustration: sensitivity-adjusted confidence intervals\n")
  print(sensitivity, digits = 6, row.names = FALSE)

  prepare_directory(directory)
  write_csv(table1, file.path(directory, "table1_empirical_results.csv"))
  write_csv(sensitivity, file.path(directory, "illustration_sensitivity.csv"))

  metadata <- list(study = "2023 PIAAC Korean empirical illustration",
    input_file = normalizePath(data_file, winslash = "/"),
    input_md5 = unname(tools::md5sum(data_file)), test_item = test_item,
    anchor_item = anchor_item, age_variable = "AGEG10LFS", reference_code = reference_code,
    focal_code = focal_code, n_reference = sum(dat$G == 0), n_focal = sum(dat$G == 1),
    plausible_values = pv_vars, weighted = FALSE, alpha = alpha, eta = eta, mu_hat = mu_hat,
    sensitivity_inputs_fixed = TRUE, confidence_intervals_clipped = FALSE,
    package_versions = vapply(c("dplyr", "sandwich"),
      function(package) as.character(utils::packageVersion(package)), character(1)),
    R_version = R.version.string)
  save_run_information(metadata, directory)

  invisible(list(table1 = table1, pv_results = pv_results, pv_differences = pv_differences,
    sensitivity = sensitivity))
}

# 9. Run the selected analysis --------------------------------------------

run_selected_tasks <- function() {
  modes <- c("simulations", "outputs", "empirical", "all", "checks", "smoke")
  if (length(RUN_MODE) != 1L || is.na(RUN_MODE) || !RUN_MODE %in% modes) {
    stop("RUN_MODE must be one of: ", paste(modes, collapse = ", "))
  }
  simulate <- RUN_MODE %in% c("simulations", "all", "smoke")
  outputs <- simulate || RUN_MODE == "outputs"
  empirical <- RUN_MODE %in% c("empirical", "all")
  replications <- if (RUN_MODE == "smoke") 20L else N_REP
  root <- if (RUN_MODE == "smoke") file.path(OUTPUT_ROOT, "smoke") else OUTPUT_ROOT
  stopifnot(length(root) == 1L, !is.na(root), nzchar(root))
  if (simulate) {
    stopifnot(length(replications) == 1L, is.finite(replications),
              replications >= 2, replications == floor(replications))
  }
  # The manuscript tables and figures use two-sided alpha = .05.
  if (ALPHA != .05) stop("This manuscript script requires ALPHA = .05.")
  if (empirical) check_empirical_inputs(EMPIRICAL_DATA_FILE)
  d1 <- file.path(root, "study1")
  d2 <- file.path(root, "study2")
  if (RUN_MODE == "outputs") {
    for (directory in c(d1, d2)) {
      path <- file.path(directory, "simulation.rds")
      if (!file.exists(path)) stop("Missing ", path, ". Run mode 'simulations' first.")
      if (!identical(readRDS(path)$metadata$alpha, .05)) {
        stop("Saved results must use alpha = .05 for the manuscript outputs.")
      }
    }
  }
  if (simulate || RUN_MODE == "checks") check_implementation()
  if (RUN_MODE == "smoke") message("20 replications per condition; not manuscript results.")
  if (simulate) simulate_study1(replications, STUDY1_SEED, d1)
  if (outputs) outputs_study1(d1)
  if (simulate) simulate_study2(replications, STUDY2_SEED, d2)
  if (outputs) outputs_study2(d2)
  if (empirical) {
    run_empirical_illustration(EMPIRICAL_DATA_FILE, file.path(root, "empirical"),
                               EMPIRICAL_ETA, ALPHA)
  }
  message("Completed: ", RUN_MODE)
  invisible(NULL)
}

run_selected_tasks()
