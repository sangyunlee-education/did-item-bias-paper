# 02_study2_falsification_diagnostic.R
# Study 2: Estimation and Inference for the Falsification Diagnostic
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
MASTER_SEED <- 20260806L + 100000L

sample_sizes <- c(500L, 1000L, 2000L)
alpha_values <- c(
  -0.10, -0.075, -0.05, -0.025, 0,
   0.025, 0.05, 0.075, 0.10
)
rho_values <- c(0, 0.3, 0.6)

make_results_dir()


# -------------------------------------------------------------------------
# Calibrate c(alpha)
# -------------------------------------------------------------------------

calibration <- data.frame(
  alpha = alpha_values,
  c_alpha = vapply(
    alpha_values,
    FUN = calibrate_c_study2,
    FUN.VALUE = numeric(1)
  )
)

write.csv(
  calibration,
  file = file.path(
    "results", "simulation", "study2_calibrated_c_alpha.csv"
  ),
  row.names = FALSE
)

cat("\nStudy 2 calibrated c(alpha) values:\n")
print(calibration, row.names = FALSE)


# -------------------------------------------------------------------------
# One Study 2 condition
# -------------------------------------------------------------------------

run_study2_condition <- function(N, alpha, rho, c_alpha,
                                 n_rep, chunk_size, seed) {
  stopifnot(N %% 2 == 0)

  set.seed(seed)

  n_group <- N / 2
  alpha_hat_all <- numeric(n_rep)
  se_alpha_all <- numeric(n_rep)

  start_indices <- seq.int(1L, n_rep, by = chunk_size)

  for (start in start_indices) {
    end <- min(start + chunk_size - 1L, n_rep)
    b <- end - start + 1L

    # Reference and focal latent-trait distributions.
    theta0 <- matrix(
      rnorm(n_group * b, mean = 0, sd = 1),
      nrow = n_group,
      ncol = b
    )

    theta1 <- matrix(
      rnorm(n_group * b, mean = -0.5, sd = 1),
      nrow = n_group,
      ncol = b
    )

    # Test item is unchanged across groups, so tau = 0.
    p_t0 <- inv_logit(theta0)
    p_t1 <- inv_logit(theta1)

    # Anchor-item shift c(alpha) is the same in both groups.
    p_a0 <- inv_logit(theta0 + c_alpha)
    p_a1 <- inv_logit(theta1 + c_alpha)

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

    alpha_hat_all[start:end] <- fit$alpha_hat
    se_alpha_all[start:end] <- fit$se_alpha
  }

  out <- summarize_mc(
    estimates = alpha_hat_all,
    ses = se_alpha_all,
    target = alpha,
    null_value = 0,
    df = N - 2
  )

  out$N <- N
  out$alpha <- alpha
  out$rho <- rho
  out$c_alpha <- c_alpha

  out[, c(
    "N", "alpha", "rho", "c_alpha",
    "Bias", "SE_SD", "Coverage", "Rejection"
  )]
}


# -------------------------------------------------------------------------
# Run all 81 conditions
# -------------------------------------------------------------------------

conditions <- expand.grid(
  N = sample_sizes,
  alpha = alpha_values,
  rho = rho_values,
  stringsAsFactors = FALSE
)

conditions$c_alpha <- vapply(
  conditions$alpha,
  FUN = function(alpha) {
    calibration$c_alpha[calibration$alpha == alpha]
  },
  FUN.VALUE = numeric(1)
)

results <- vector("list", nrow(conditions))

cat("\nRunning Study 2:", nrow(conditions), "conditions x",
    N_REP, "replications\n\n")

for (i in seq_len(nrow(conditions))) {
  cond <- conditions[i, ]

  cat(
    sprintf(
      "[Study 2 %02d/%02d] N=%d, alpha=%.3f, rho=%.1f\n",
      i, nrow(conditions), cond$N, cond$alpha, cond$rho
    )
  )

  results[[i]] <- run_study2_condition(
    N = cond$N,
    alpha = cond$alpha,
    rho = cond$rho,
    c_alpha = cond$c_alpha,
    n_rep = N_REP,
    chunk_size = CHUNK_SIZE,
    seed = MASTER_SEED + i
  )
}

study2_results <- do.call(rbind, results)

write.csv(
  study2_results,
  file = file.path(
    "results", "simulation", "study2_condition_results.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    "results", "simulation", "sessionInfo_study2.txt"
  )
)

cat("\nStudy 2 complete.\n")
cat("Saved: results/simulation/study2_condition_results.csv\n")
