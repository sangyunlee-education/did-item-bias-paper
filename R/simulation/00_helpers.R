# 00_helpers.R
# Helper functions for the Monte Carlo simulation studies
# "Identifying Item Bias Without Conditioning: A Difference-in-Differences Approach"
#
# This file uses base R only.

inv_logit <- function(x) {
  plogis(x)
}


# -------------------------------------------------------------------------
# Latent-trait distributions
# -------------------------------------------------------------------------

# Parameters for W ~ SN(0, 1, lambda), with lambda = 2.17.
# A skew-normal draw can be generated as
#   W = delta * |Z0| + sqrt(1 - delta^2) * Z1,
# where Z0 and Z1 are independent standard normal variables.
skew_lambda <- 2.17
skew_delta <- skew_lambda / sqrt(1 + skew_lambda^2)
skew_mean_w <- skew_delta * sqrt(2 / pi)
skew_sd_w <- sqrt(1 - 2 * skew_delta^2 / pi)

rtheta_focal <- function(n, distribution) {
  if (distribution == "Normal") {
    return(rnorm(n, mean = -0.5, sd = 1))
  }

  if (distribution == "Skewed") {
    z0 <- rnorm(n)
    z1 <- rnorm(n)
    w <- skew_delta * abs(z0) + sqrt(1 - skew_delta^2) * z1

    # Standardize W, then set the focal-group mean to -0.5 and SD to 1.
    return(-0.5 + (w - skew_mean_w) / skew_sd_w)
  }

  if (distribution == "Bimodal") {
    component <- rbinom(n, size = 1, prob = 0.5)
    means <- ifelse(component == 0, -1.3, 0.3)
    return(rnorm(n, mean = means, sd = 0.6))
  }

  stop("Unknown focal-group distribution: ", distribution)
}


# -------------------------------------------------------------------------
# Numerical calibration
# -------------------------------------------------------------------------

# Study 1:
# Find delta_T such that
# E_focal[P_T(1, theta) - P_T(0, theta)] = tau,
# where P_T(0, theta) = logistic(theta) and
#       P_T(1, theta) = logistic(theta + delta_T).
calibrate_delta_study1 <- function(tau, distribution) {
  if (tau == 0) {
    return(0)
  }

  target_difference <- function(delta_T) {
    if (distribution == "Normal") {
      value <- integrate(
        f = function(theta) {
          (inv_logit(theta + delta_T) - inv_logit(theta)) *
            dnorm(theta, mean = -0.5, sd = 1)
        },
        lower = -Inf,
        upper = Inf,
        rel.tol = 1e-11
      )$value
    } else if (distribution == "Skewed") {
      # Integrate over the original skew-normal variable W.
      value <- integrate(
        f = function(w) {
          theta <- -0.5 + (w - skew_mean_w) / skew_sd_w
          d_w <- 2 * dnorm(w) * pnorm(skew_lambda * w)

          (inv_logit(theta + delta_T) - inv_logit(theta)) * d_w
        },
        lower = -Inf,
        upper = Inf,
        rel.tol = 1e-11,
        subdivisions = 1000L
      )$value
    } else if (distribution == "Bimodal") {
      value <- integrate(
        f = function(theta) {
          d_theta <- 0.5 * dnorm(theta, mean = -1.3, sd = 0.6) +
            0.5 * dnorm(theta, mean = 0.3, sd = 0.6)

          (inv_logit(theta + delta_T) - inv_logit(theta)) * d_theta
        },
        lower = -Inf,
        upper = Inf,
        rel.tol = 1e-11
      )$value
    } else {
      stop("Unknown focal-group distribution: ", distribution)
    }

    value - tau
  }

  uniroot(
    f = target_difference,
    interval = c(-3, 3),
    tol = 1e-11
  )$root
}


# Study 2:
# Find c(alpha) such that
# E[Y_T - Y_A | G = 0] = alpha,
# with
#   P_T(theta) = logistic(theta)
#   P_A(theta) = logistic(theta + c(alpha)).
calibrate_c_study2 <- function(alpha) {
  if (alpha == 0) {
    return(0)
  }

  target_difference <- function(c_alpha) {
    value <- integrate(
      f = function(theta) {
        (inv_logit(theta) - inv_logit(theta + c_alpha)) * dnorm(theta)
      },
      lower = -Inf,
      upper = Inf,
      rel.tol = 1e-11
    )$value

    value - alpha
  }

  uniroot(
    f = target_difference,
    interval = c(-3, 3),
    tol = 1e-11
  )$root
}


# -------------------------------------------------------------------------
# HC3 inference for E(D | G) = alpha + beta G, where D = Y_T - Y_A
# -------------------------------------------------------------------------

# With an intercept and a binary group indicator, the regression is
# equivalent to comparing the two group means. The leverage of each
# observation in group g is 1 / n_g. The formulas below are exactly the
# HC3 formulas for this two-group regression and avoid repeatedly fitting
# hundreds of thousands of lm() models.
#
# d0 and d1 are matrices:
#   rows    = respondents within a group
#   columns = Monte Carlo replications in the current chunk.
hc3_two_group <- function(d0, d1) {
  n0 <- nrow(d0)
  n1 <- nrow(d1)

  mean0 <- colMeans(d0)
  mean1 <- colMeans(d1)

  # Sum of squared residuals within each group.
  sse0 <- colSums(d0^2) - n0 * mean0^2
  sse1 <- colSums(d1^2) - n1 * mean1^2

  # Guard against tiny negative values from floating-point arithmetic.
  sse0 <- pmax(sse0, 0)
  sse1 <- pmax(sse1, 0)

  # HC3 variance of each group mean.
  v0 <- sse0 / (n0 - 1)^2
  v1 <- sse1 / (n1 - 1)^2

  alpha_hat <- mean0
  beta_hat <- mean1 - mean0

  se_alpha <- sqrt(v0)
  se_beta <- sqrt(v0 + v1)

  list(
    alpha_hat = alpha_hat,
    beta_hat = beta_hat,
    se_alpha = se_alpha,
    se_beta = se_beta
  )
}


# -------------------------------------------------------------------------
# Condition-level summaries
# -------------------------------------------------------------------------

summarize_mc <- function(estimates, ses, target, null_value, df) {
  critical <- qt(0.975, df = df)

  lower <- estimates - critical * ses
  upper <- estimates + critical * ses

  # When a standard error is numerically zero, treat the test statistic
  # as zero if the estimate also equals the null and infinite otherwise.
  test_stat <- ifelse(
    ses > 0,
    (estimates - null_value) / ses,
    ifelse(estimates == null_value, 0, Inf)
  )

  data.frame(
    Bias = mean(estimates) - target,
    SE_SD = mean(ses) / sd(estimates),
    Coverage = mean(lower <= target & upper >= target),
    Rejection = mean(abs(test_stat) > critical)
  )
}


# -------------------------------------------------------------------------
# Utility
# -------------------------------------------------------------------------

make_results_dir <- function() {
  dir.create(
    file.path("results", "simulation"),
    recursive = TRUE,
    showWarnings = FALSE
  )
}
