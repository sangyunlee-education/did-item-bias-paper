# 01_study1_exact_equivalence.R
# Study 1: Finite-Sample Inference Under Exact Identification
#
# Run this script from the repository root.
# It sources R/simulation/00_helpers.R and writes results to
# results/simulation/.
#
# Full paper run: N_REP <- 5000L
# For a quick code check, temporarily use N_REP <- 50L or 100L.

source(file.path("R", "simulation", "00_helpers.R"))

N_REP <- 5000L
CHUNK_SIZE <- 100L
MASTER_SEED <- 20260806L

sample_sizes <- c(500L, 1000L, 2000L)
distributions <- c("Normal", "Skewed", "Bimodal")
tau_values <- c(0, -0.05, -0.10)
rho_values <- c(0, 0.3, 0.6)

make_results_dir()


# -------------------------------------------------------------------------
# Calibrate delta_T
# -------------------------------------------------------------------------

calibration <- expand.grid(
  Distribution = distributions,
  tau = tau_values,
  stringsAsFactors = FALSE
)

calibration$delta_T <- mapply(
  FUN = calibrate_delta_study1,
  tau = calibration$tau,
  distribution = calibration$Distribution
)

write.csv(
  calibration,
  file = file.path(
    "results", "simulation", "study1_calibrated_delta_T.csv"
  ),
  row.names = FALSE
)

cat("\nStudy 1 calibrated delta_T values:\n")
print(calibration, row.names = FALSE)


# -------------------------------------------------------------------------
# One Study 1 condition
# -------------------------------------------------------------------------

run_study1_condition <- function(N, distribution, tau, rho,
                                 delta_T, n_rep, chunk_size, seed) {
  stopifnot(N %% 2 == 0)

  set.seed(seed)

  n_group <- N / 2
  beta_hat_all <- numeric(n_rep)
  se_beta_all <- numeric(n_rep)

  start_indices <- seq.int(1L, n_rep, by = chunk_size)

  for (start in start_indices) {
    end <- min(start + chunk_size - 1L, n_rep)
    b <- end - start + 1L

    # Reference group: Theta ~ N(0, 1)
    theta0 <- matrix(
      rnorm(n_group * b, mean = 0, sd = 1),
      nrow = n_group,
      ncol = b
    )

    # Focal group: distribution specified by the condition
    theta1 <- matrix(
      rtheta_focal(n_group * b, distribution = distribution),
      nrow = n_group,
      ncol = b
    )

    # Exact item response function equivalence in the reference condition:
    # P_A(theta) = P_T(0, theta) = logistic(theta)
    p_a0 <- inv_logit(theta0)
    p_t0 <- inv_logit(theta0)

    p_a1 <- inv_logit(theta1)
    p_t1 <- inv_logit(theta1 + delta_T)

    # Gaussian-copula construction for residual dependence.
    u_a0 <- matrix(rnorm(n_group * b), nrow = n_group, ncol = b)
    eps0 <- matrix(rnorm(n_group * b), nrow = n_group, ncol = b)
    u_t0 <- rho * u_a0 + sqrt(1 - rho^2) * eps0

    u_a1 <- matrix(rnorm(n_group * b), nrow = n_group, ncol = b)
    eps1 <- matrix(rnorm(n_group * b), nrow = n_group, ncol = b)
    u_t1 <- rho * u_a1 + sqrt(1 - rho^2) * eps1

    y_a0 <- (pnorm(u_a0) <= p_a0)
    y_t0 <- (pnorm(u_t0) <= p_t0)

    y_a1 <- (pnorm(u_a1) <= p_a1)
    y_t1 <- (pnorm(u_t1) <= p_t1)

    d0 <- y_t0 - y_a0
    d1 <- y_t1 - y_a1

    fit <- hc3_two_group(d0, d1)

    beta_hat_all[start:end] <- fit$beta_hat
    se_beta_all[start:end] <- fit$se_beta
  }

  out <- summarize_mc(
    estimates = beta_hat_all,
    ses = se_beta_all,
    target = tau,
    null_value = 0,
    df = N - 2
  )

  out$N <- N
  out$Distribution <- distribution
  out$tau <- tau
  out$rho <- rho
  out$delta_T <- delta_T

  out[, c(
    "N", "Distribution", "tau", "rho", "delta_T",
    "Bias", "SE_SD", "Coverage", "Rejection"
  )]
}


# -------------------------------------------------------------------------
# Run all 81 conditions
# -------------------------------------------------------------------------

conditions <- expand.grid(
  N = sample_sizes,
  Distribution = distributions,
  tau = tau_values,
  rho = rho_values,
  stringsAsFactors = FALSE
)

conditions$delta_T <- mapply(
  FUN = function(distribution, tau) {
    calibration$delta_T[
      calibration$Distribution == distribution &
        calibration$tau == tau
    ]
  },
  distribution = conditions$Distribution,
  tau = conditions$tau
)

results <- vector("list", nrow(conditions))

cat("\nRunning Study 1:", nrow(conditions), "conditions x",
    N_REP, "replications\n\n")

for (i in seq_len(nrow(conditions))) {
  cond <- conditions[i, ]

  cat(
    sprintf(
      "[Study 1 %02d/%02d] N=%d, Distribution=%s, tau=%.3f, rho=%.1f\n",
      i, nrow(conditions), cond$N, cond$Distribution, cond$tau, cond$rho
    )
  )

  results[[i]] <- run_study1_condition(
    N = cond$N,
    distribution = cond$Distribution,
    tau = cond$tau,
    rho = cond$rho,
    delta_T = cond$delta_T,
    n_rep = N_REP,
    chunk_size = CHUNK_SIZE,
    seed = MASTER_SEED + i
  )
}

study1_results <- do.call(rbind, results)

write.csv(
  study1_results,
  file = file.path(
    "results", "simulation", "study1_condition_results.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    "results", "simulation", "sessionInfo_study1.txt"
  )
)

cat("\nStudy 1 complete.\n")
cat("Saved: results/simulation/study1_condition_results.csv\n")
