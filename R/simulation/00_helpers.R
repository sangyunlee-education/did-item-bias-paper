# 00_helpers.R
# Helper functions for the Monte Carlo simulation studies
# "Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach"
#
# IMPORTANT FOR REPRODUCTION
# This file preserves the random-number generation logic used in the original
# simulation script that generated the results reported in the manuscript.
# In particular, it preserves the original base seed, block size, latent-trait
# generators, binary-response generator, and HC3 calculations.

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

SETTINGS <- list(
  seed = 20260806L,
  n_rep = 5000L,
  quick_test = FALSE,
  n_cores = max(1L, parallel::detectCores(logical = TRUE) - 1L),
  chunk_size = 100L,
  output_dir = file.path("results", "simulation")
)

# For a quick code check only, change quick_test above to TRUE.
# The manuscript results use quick_test = FALSE and n_rep = 5000.
if (isTRUE(SETTINGS$quick_test)) {
  SETTINGS$n_rep <- 100L
  SETTINGS$n_cores <- min(2L, SETTINGS$n_cores)
  SETTINGS$chunk_size <- 25L
}

dir.create(SETTINGS$output_dir, recursive = TRUE, showWarnings = FALSE)

# Baseline item parameters
A_T <- 1
B_T <- 0

# Study 1 factors
N_LEVELS <- c(500L, 1000L, 2000L)
TAU_LEVELS <- c(0, -0.05, -0.10)
RHO_LEVELS <- c(0, 0.3, 0.6)
DIST_LEVELS <- c("Normal", "Skewed", "Bimodal")

# Study 2 factors
ALPHA_LEVELS <- c(-0.10, -0.075, -0.05, -0.025, 0,
                  0.025, 0.05, 0.075, 0.10)
MU1_DIAGNOSTIC <- -0.50

# Skew-normal constants used in Study 1
SN_SHAPE <- 2.17
SN_DELTA <- SN_SHAPE / sqrt(1 + SN_SHAPE^2)
SN_MEAN <- SN_DELTA * sqrt(2 / pi)
SN_SD <- sqrt(1 - 2 * SN_DELTA^2 / pi)

# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------

integrate_safe <- function(fun) {
  integrate(
    fun,
    lower = -Inf,
    upper = Inf,
    subdivisions = 2000L,
    rel.tol = 1e-10,
    abs.tol = 1e-10,
    stop.on.error = TRUE
  )$value
}

format_decimal <- function(x, digits = 3L, drop_leading_zero = TRUE) {
  out <- formatC(x, format = "f", digits = digits)
  out[abs(x) < 0.5 * 10^(-digits)] <- "0"
  if (drop_leading_zero) {
    out <- sub("^0\\.", ".", out)
    out <- sub("^-0\\.", "-.", out)
  }
  out
}

format_range <- function(low, high, digits = 3L) {
  paste0(
    format_decimal(low, digits = digits),
    "--",
    format_decimal(high, digits = digits)
  )
}

parallel_lapply_conditions <- function(condition_list, fun) {
  if (SETTINGS$n_cores <= 1L || length(condition_list) <= 1L) {
    return(lapply(condition_list, fun))
  }

  if (.Platform$OS.type == "windows") {
    cluster <- parallel::makeCluster(SETTINGS$n_cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterExport(
      cluster,
      varlist = ls(envir = .GlobalEnv),
      envir = .GlobalEnv
    )
    return(parallel::parLapply(cluster, condition_list, fun))
  }

  parallel::mclapply(
    condition_list,
    fun,
    mc.cores = SETTINGS$n_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
}

# -----------------------------------------------------------------------------
# Latent-trait distributions
# -----------------------------------------------------------------------------

r_theta_reference <- function(n) {
  rnorm(n, mean = 0, sd = 1)
}

r_theta_focal_study1 <- function(n, distribution) {
  if (distribution == "Normal") {
    return(rnorm(n, mean = -0.5, sd = 1))
  }

  if (distribution == "Skewed") {
    u0 <- rnorm(n)
    u1 <- rnorm(n)
    z <- SN_DELTA * abs(u0) + sqrt(1 - SN_DELTA^2) * u1
    return(-0.5 + (z - SN_MEAN) / SN_SD)
  }

  if (distribution == "Bimodal") {
    component <- rbinom(n, size = 1, prob = 0.5)
    means <- ifelse(component == 1, -1.3, 0.3)
    return(rnorm(n, mean = means, sd = 0.6))
  }

  stop("Unknown Study 1 focal distribution: ", distribution)
}

r_theta_focal_study2 <- function(n, mu1) {
  rnorm(n, mean = mu1, sd = 1)
}

f0_density <- function(theta) {
  dnorm(theta, mean = 0, sd = 1)
}

f1_density_study1 <- function(theta, distribution) {
  if (distribution == "Normal") {
    return(dnorm(theta, mean = -0.5, sd = 1))
  }

  if (distribution == "Skewed") {
    # Theta = -0.5 + (Z - E[Z]) / SD[Z], Z ~ SN(0, 1, SN_SHAPE)
    z <- SN_MEAN + SN_SD * (theta + 0.5)
    f_z <- 2 * dnorm(z) * pnorm(SN_SHAPE * z)
    return(SN_SD * f_z)
  }

  if (distribution == "Bimodal") {
    return(
      0.5 * dnorm(theta, mean = -1.3, sd = 0.6) +
        0.5 * dnorm(theta, mean = 0.3, sd = 0.6)
    )
  }

  stop("Unknown Study 1 focal distribution: ", distribution)
}

# -----------------------------------------------------------------------------
# Item response functions and population targets
# -----------------------------------------------------------------------------

p_test_reference <- function(theta) {
  plogis(A_T * (theta - B_T))
}

p_test_focal <- function(theta, delta_t) {
  plogis(A_T * (theta - B_T) + delta_t)
}

p_anchor_shift <- function(theta, shift_c) {
  plogis(A_T * (theta - B_T) + shift_c)
}

find_delta_t <- function(target_tau, distribution) {
  if (abs(target_tau) < 1e-14) {
    return(0)
  }

  objective <- function(delta_t) {
    integrate_safe(function(theta) {
      (p_test_focal(theta, delta_t) - p_test_reference(theta)) *
        f1_density_study1(theta, distribution)
    }) - target_tau
  }

  uniroot(objective, interval = c(-15, 15), tol = 1e-11)$root
}

alpha_from_shift <- function(shift_c) {
  integrate_safe(function(theta) {
    (p_test_reference(theta) - p_anchor_shift(theta, shift_c)) *
      f0_density(theta)
  })
}

find_anchor_shift <- function(target_alpha) {
  if (abs(target_alpha) < 1e-14) {
    return(0)
  }

  objective <- function(shift_c) {
    alpha_from_shift(shift_c) - target_alpha
  }

  uniroot(objective, interval = c(-15, 15), tol = 1e-11)$root
}

# -----------------------------------------------------------------------------
# Binary item generation and HC3 inference
# -----------------------------------------------------------------------------

generate_binary_pair <- function(p_a, p_t, rho) {
  if (length(p_a) != length(p_t)) {
    stop("p_a and p_t must have the same length.")
  }
  if (abs(rho) >= 1) {
    stop("rho must lie strictly between -1 and 1.")
  }

  # Do not change the order of these random draws if exact reproduction of
  # the manuscript simulation is desired.
  z_a <- rnorm(length(p_a))
  z_t <- rho * z_a + sqrt(1 - rho^2) * rnorm(length(p_a))

  list(
    y_a = as.integer(pnorm(z_a) <= p_a),
    y_t = as.integer(pnorm(z_t) <= p_t)
  )
}

row_block_statistics <- function(values, group_size, n_rep_block) {
  value_matrix <- matrix(values, nrow = group_size, ncol = n_rep_block)
  means <- colMeans(value_matrix)
  sum_squares <- colSums(value_matrix^2)
  variances <- pmax(
    (sum_squares - group_size * means^2) / (group_size - 1),
    0
  )
  list(mean = means, variance = variances)
}

simulate_blocks <- function(
    n0,
    n1,
    rho,
    r_theta_focal,
    p_a_reference,
    p_t_reference,
    p_a_focal,
    p_t_focal) {

  n_rep <- SETTINGS$n_rep
  alpha_hat <- numeric(n_rep)
  beta_hat <- numeric(n_rep)
  se_alpha <- numeric(n_rep)
  se_beta <- numeric(n_rep)
  p_alpha <- numeric(n_rep)
  p_beta <- numeric(n_rep)
  ci_alpha_low <- numeric(n_rep)
  ci_alpha_high <- numeric(n_rep)
  ci_beta_low <- numeric(n_rep)
  ci_beta_high <- numeric(n_rep)

  residual_df <- n0 + n1 - 2L
  critical_value <- qt(0.975, df = residual_df)
  starts <- seq.int(1L, n_rep, by = SETTINGS$chunk_size)

  for (start in starts) {
    end <- min(start + SETTINGS$chunk_size - 1L, n_rep)
    indices <- start:end
    block_size <- length(indices)

    # Do not change this draw order if exact reproduction is desired.
    theta0 <- r_theta_reference(n0 * block_size)
    theta1 <- r_theta_focal(n1 * block_size)

    pair0 <- generate_binary_pair(
      p_a = p_a_reference(theta0),
      p_t = p_t_reference(theta0),
      rho = rho
    )
    pair1 <- generate_binary_pair(
      p_a = p_a_focal(theta1),
      p_t = p_t_focal(theta1),
      rho = rho
    )

    d0 <- row_block_statistics(
      pair0$y_t - pair0$y_a,
      group_size = n0,
      n_rep_block = block_size
    )
    d1 <- row_block_statistics(
      pair1$y_t - pair1$y_a,
      group_size = n1,
      n_rep_block = block_size
    )

    alpha_block <- d0$mean
    beta_block <- d1$mean - d0$mean

    # In the saturated two-group regression, the HC3 variance of a group
    # mean is s_g^2/(n_g - 1). The slope variance is the sum of the two
    # group-mean variances.
    variance0_hc3 <- d0$variance / (n0 - 1)
    variance1_hc3 <- d1$variance / (n1 - 1)
    se_alpha_block <- sqrt(variance0_hc3)
    se_beta_block <- sqrt(variance0_hc3 + variance1_hc3)

    t_alpha <- ifelse(
      se_alpha_block > 0,
      alpha_block / se_alpha_block,
      ifelse(abs(alpha_block) < 1e-14, 0, sign(alpha_block) * Inf)
    )
    t_beta <- ifelse(
      se_beta_block > 0,
      beta_block / se_beta_block,
      ifelse(abs(beta_block) < 1e-14, 0, sign(beta_block) * Inf)
    )

    alpha_hat[indices] <- alpha_block
    beta_hat[indices] <- beta_block
    se_alpha[indices] <- se_alpha_block
    se_beta[indices] <- se_beta_block
    p_alpha[indices] <- 2 * pt(-abs(t_alpha), df = residual_df)
    p_beta[indices] <- 2 * pt(-abs(t_beta), df = residual_df)
    ci_alpha_low[indices] <- alpha_block - critical_value * se_alpha_block
    ci_alpha_high[indices] <- alpha_block + critical_value * se_alpha_block
    ci_beta_low[indices] <- beta_block - critical_value * se_beta_block
    ci_beta_high[indices] <- beta_block + critical_value * se_beta_block
  }

  list(
    alpha_hat = alpha_hat,
    beta_hat = beta_hat,
    se_alpha = se_alpha,
    se_beta = se_beta,
    p_alpha = p_alpha,
    p_beta = p_beta,
    ci_alpha_low = ci_alpha_low,
    ci_alpha_high = ci_alpha_high,
    ci_beta_low = ci_beta_low,
    ci_beta_high = ci_beta_high
  )
}

# -----------------------------------------------------------------------------
# Study-specific condition functions
# -----------------------------------------------------------------------------

simulate_study1_condition <- function(condition) {
  set.seed(as.integer(condition$condition_seed))

  n_total <- as.integer(condition$N)
  n0 <- n_total %/% 2L
  n1 <- n_total - n0
  distribution <- as.character(condition$distribution)
  rho <- as.numeric(condition$rho)
  tau <- as.numeric(condition$tau)
  delta_t <- as.numeric(condition$delta_t)

  simulation <- simulate_blocks(
    n0 = n0,
    n1 = n1,
    rho = rho,
    r_theta_focal = function(n) r_theta_focal_study1(n, distribution),
    p_a_reference = p_test_reference,
    p_t_reference = p_test_reference,
    p_a_focal = p_test_reference,
    p_t_focal = function(theta) p_test_focal(theta, delta_t)
  )

  empirical_sd <- sd(simulation$beta_hat)

  data.frame(
    N = n_total,
    rho = rho,
    distribution = distribution,
    tau = tau,
    delta_t = delta_t,
    replications = SETTINGS$n_rep,
    mc_mean_did = mean(simulation$beta_hat),
    mc_bias = mean(simulation$beta_hat) - tau,
    rmse = sqrt(mean((simulation$beta_hat - tau)^2)),
    mean_se = mean(simulation$se_beta),
    empirical_sd = empirical_sd,
    se_sd = mean(simulation$se_beta) / empirical_sd,
    coverage = mean(
      simulation$ci_beta_low <= tau & simulation$ci_beta_high >= tau
    ),
    reject_beta = mean(simulation$p_beta < 0.05),
    mean_alpha = mean(simulation$alpha_hat),
    mean_se_alpha = mean(simulation$se_alpha),
    empirical_sd_alpha = sd(simulation$alpha_hat),
    reject_alpha = mean(simulation$p_alpha < 0.05),
    coverage_alpha = mean(
      simulation$ci_alpha_low <= 0 & simulation$ci_alpha_high >= 0
    ),
    stringsAsFactors = FALSE
  )
}

simulate_study2_condition <- function(condition) {
  set.seed(as.integer(condition$condition_seed))

  n_total <- as.integer(condition$N)
  n0 <- n_total %/% 2L
  n1 <- n_total - n0
  rho <- as.numeric(condition$rho)
  mu1 <- as.numeric(condition$mu1)
  target_alpha <- as.numeric(condition$alpha)
  shift_c <- as.numeric(condition$shift_c)

  anchor_function <- function(theta) p_anchor_shift(theta, shift_c)

  simulation <- simulate_blocks(
    n0 = n0,
    n1 = n1,
    rho = rho,
    r_theta_focal = function(n) r_theta_focal_study2(n, mu1),
    p_a_reference = anchor_function,
    p_t_reference = p_test_reference,
    p_a_focal = anchor_function,
    p_t_focal = p_test_reference
  )

  empirical_sd_alpha <- sd(simulation$alpha_hat)

  data.frame(
    N = n_total,
    rho = rho,
    mu1 = mu1,
    alpha = target_alpha,
    shift_c = shift_c,
    replications = SETTINGS$n_rep,
    mc_mean_alpha = mean(simulation$alpha_hat),
    mc_bias_alpha = mean(simulation$alpha_hat) - target_alpha,
    mean_se_alpha = mean(simulation$se_alpha),
    empirical_sd_alpha = empirical_sd_alpha,
    se_sd_alpha = mean(simulation$se_alpha) / empirical_sd_alpha,
    coverage_alpha = mean(
      simulation$ci_alpha_low <= target_alpha &
        simulation$ci_alpha_high >= target_alpha
    ),
    reject_alpha = mean(simulation$p_alpha < 0.05),
    stringsAsFactors = FALSE
  )
}
