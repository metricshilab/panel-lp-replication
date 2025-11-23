# =============================================================================
# 02-02-regression-methods.R
# Purpose: Implement regression methods for panel local projection estimation
# =============================================================================

# =============================================================================
# OLS Regression
# =============================================================================

# Perform OLS regression with dynamic formula construction
ols_regression <- function(data, dep, indep, dep_lag = NULL, lag = FALSE, indep2 = NULL, vars_num = 2) {
  formula_parts <- c(indep, if (lag) dep_lag, if (vars_num == 3) indep2)
  formula <- as.formula(paste(dep, "~", paste(formula_parts, collapse = " + "), "- 1"))
  model <- lm(formula, data = data)
  return(as.numeric(coef(model)))
}

# =============================================================================
# Estimation Methods
# =============================================================================

# Compute Split-Panel Jackknife estimates
compute_spj_estimates <- function(data, midpoint, lag = FALSE, vars_num = 2) {
  # Full-sample estimate
  full_beta <- ols_regression(data, "dep_lead_h_demean", "indep_demean", "dep_lag_demean", lag, "indep2_demean", vars_num)

  # Subsample estimates for periods A and B
  beta_a <- ols_regression(data[time <= midpoint], "dep_lead_h_demean_ab", "indep_demean_ab", "dep_lag_demean_ab", lag, "indep2_demean_ab", vars_num)
  beta_b <- ols_regression(data[time > midpoint], "dep_lead_h_demean_ab", "indep_demean_ab", "dep_lag_demean_ab", lag, "indep2_demean_ab", vars_num)

  # SPJ bias-corrected estimate
  return(2 * full_beta - 0.5 * (beta_a + beta_b))
}

# Compute Fixed Effects estimates
compute_fe_estimates <- function(data, lag = FALSE, vars_num = 2) {
  return(ols_regression(data, "dep_lead_h_demean", "indep_demean", "dep_lag_demean", lag, "indep2_demean", vars_num))
}

# Compute Difference GMM estimates
compute_dgmm_estimates <- function(data) {
  model <- AER::ivreg(
    dep_lead_h_diff ~ indep_diff + dep_lag_diff - 1 |
      iv_indep + iv_dep - 1,
    data = data
  )
  return(list(coefficients = coef(model), residuals = resid(model)))
}
