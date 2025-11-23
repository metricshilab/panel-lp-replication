# =============================================================================
# 02-estimate.R
# Purpose: Main estimation functions for panel local projection methods
# =============================================================================

# Load dependencies
source("scripts/02-01-data-preparation.R")
source("scripts/02-02-regression-methods.R")
source("scripts/02-03-error-calculation.R")

# =============================================================================
# DGMM Estimation
# =============================================================================

# Estimate IRF using Difference GMM
estimate_dgmm <- function(data, num_periods, num_individuals,
                          horizon = 10, time_effects = TRUE, lag = FALSE,
                          cluster = FALSE, dif = FALSE, vars_num = 2,
                          cluster_method = "driscoll_kraay") {
  results <- lapply(0:horizon, function(h) {
    # Prepare and difference data
    data_horizon <- prepare_data(data, h, lag = lag, dif = dif)
    T_h <- num_periods - h
    data_diff <- difference_data(data_horizon)

    # Compute estimates
    estimate_results <- compute_dgmm_estimates(data_diff)
    beta <- as.numeric(estimate_results$coefficients[1])
    residuals <- estimate_results$residuals

    # Compute standard errors
    data_longer <- longer_data_se_dgmm(data_diff, residuals)
    se <- calculate_se_dgmm(data_longer, N = num_individuals, cluster)

    return(list(Horizon = h, beta = beta, se = se, estimator = "DGMM"))
  })

  return(rbindlist(results))
}

# =============================================================================
# FE Estimation
# =============================================================================

# Estimate IRF using Fixed Effects
estimate_fe <- function(data, num_periods, num_individuals,
                        horizon = 10, time_effects = TRUE, lag = FALSE,
                        cluster = FALSE, dif = FALSE, vars_num = 2,
                        cluster_method = "driscoll_kraay") {
  # Determine horizon range
  if (vars_num == 3 | dif == TRUE) {
    h_range <- 1:horizon
  } else {
    h_range <- 0:horizon
  }

  results <- lapply(h_range, function(h) {
    # Prepare and demean data
    data_horizon <- prepare_data(data, h, lag = lag, dif = dif)
    data_demean <- demean_data(data_horizon, lag = lag, vars_num = vars_num)
    rm(data_horizon)
    gc()

    # Compute estimates
    beta <- compute_fe_estimates(data_demean, lag = lag, vars_num = vars_num)

    if (vars_num == 2) {
      # Compute standard errors
      data_longer <- longer_data_se_fe(data_demean, beta, lag = lag)
      rm(data_demean)
      gc()

      se <- calculate_se_fe(data_longer, N = num_individuals, lag = lag, cluster = cluster,
                            cluster_method = cluster_method)
      rm(data_longer)
      gc()

      if (lag) {
        return(list(Horizon = h, beta = beta[[1]], se = se, estimator = "FE"))
      } else {
        return(list(Horizon = h, beta = beta, se = se, estimator = "FE"))
      }
    } else {
      return(list(Horizon = h, beta = beta[[1]], estimator = "FE"))
    }
  })

  return(rbindlist(results))
}

# =============================================================================
# SPJ Estimation
# =============================================================================

# Estimate IRF using Split-Panel Jackknife
estimate_spj <- function(data, num_periods, num_individuals,
                         horizon = 10, time_effects = FALSE, lag = FALSE,
                         cluster = FALSE, dif = FALSE, vars_num = 2,
                         cluster_method = "driscoll_kraay") {
  # Determine horizon range
  if (vars_num == 3 | dif == TRUE) {
    h_range <- 1:horizon
  } else {
    h_range <- 0:horizon
  }

  results <- lapply(h_range, function(h) {
    # Prepare data and determine midpoint
    data_horizon <- prepare_data(data, h, lag = lag, dif = dif)
    T_h <- uniqueN(data_horizon$time)
    midpoint <- unique(data_horizon$time)[ceiling(T_h / 2)]

    # Demean data and compute estimates
    data_demean <- demean_data(data_horizon, midpoint, lag = lag, vars_num = vars_num)
    rm(data_horizon)
    gc()

    beta <- compute_spj_estimates(data_demean, midpoint, lag = lag, vars_num = vars_num)

    if (vars_num == 2) {
      # Compute standard errors
      data_longer <- longer_data_se_spj(data_demean, midpoint, beta, lag = lag)
      rm(data_demean)
      gc()

      se <- calculate_se_spj(data_longer, num_individuals, midpoint, beta, lag = lag,
                             cluster = cluster, cluster_method = cluster_method)
      rm(data_longer)
      gc()

      if (lag) {
        return(list(Horizon = h, beta = beta[[1]], se = se, estimator = "SPJ"))
      } else {
        return(list(Horizon = h, beta = beta, se = se, estimator = "SPJ"))
      }
    } else {
      return(list(Horizon = h, beta = beta[[1]], estimator = "SPJ"))
    }
  })

  return(rbindlist(results))
}

# =============================================================================
# DB Estimation (De-biased)
# =============================================================================

# Estimate IRF using De-biased estimator
estimate_db <- function(data, num_periods, num_individuals,
                        horizon = 10, time_effects = TRUE, lag = FALSE,
                        cluster = FALSE, dif = FALSE, vars_num = 2,
                        cluster_method = "driscoll_kraay") {
  # Get FE estimates as base
  results_fe <- estimate_fe(
    data, num_periods, num_individuals,
    horizon = horizon, time_effects = time_effects, lag = lag,
    cluster = cluster, dif = dif, vars_num = vars_num)
  setDT(results_fe)
  results_db <- results_fe[, .(Horizon, se, estimator)]
  results_db[, estimator := rep("DB", times = (1 + horizon))]

  # Estimate AR coefficient
  T0 = num_periods
  x.f <- data[time >= 2]$indep
  x.l <- data[time < T0]$indep
  AR.fit <- summary(lm(x.f~x.l))
  rho.hat <- AR.fit$coefficients[2]
  s.u <- sqrt(mean(AR.fit$residuals^2))

  # Get initial coefficient estimate
  h = 0
  data_horizon <- prepare_data(data, h, lag = lag, dif = dif)
  data_demean <- demean_data(data_horizon, lag = lag, vars_num = vars_num)
  rm(data_horizon)
  gc()

  y.vec <- data_demean$dep_lead_h_demean
  x.vec <- data_demean$indep_demean
  all.fit <- summary(lm(y.vec~x.vec-1))
  b0.hat <- all.fit$coefficients[1]

  # Apply bias correction
  f.rho <- s.u^2*( (1-rho.hat^(0:horizon)) - (0:horizon/(T0-0:horizon)))/(1-rho.hat)^2
  sx = sqrt(mean(x.vec^2))

  results_db[, beta := results_fe$beta + b0.hat * f.rho / ( (T0-0:horizon) * sx^2)]
  return(results_db)
}

# =============================================================================
# LE Estimation (Leave-out)
# =============================================================================

# Estimate IRF using Leave-out estimator
estimate_le <- function(data, num_periods, num_individuals,
                        horizon = 10, time_effects = TRUE, lag = FALSE,
                        cluster = FALSE, dif = FALSE, vars_num = 2,
                        cluster_method = "driscoll_kraay") {
  # Get FE estimates as base
  results_fe <- estimate_fe(
    data, num_periods, num_individuals,
    horizon = horizon, time_effects = time_effects, lag = lag,
    cluster = cluster, dif = dif, vars_num = vars_num)
  setDT(results_fe)

  results_le <- results_fe[, .(Horizon, se, estimator)]
  results_le[, estimator := rep("LE", times = (1 + horizon))]
  T0 = num_periods

  # Estimate rho using Newton's method
  data_func <- copy(data)
  data_func[, indep_lag := shift(indep, 1, type = "lag"), by = id]
  data_func <- na.omit(data_func, cols = "indep_lag")

  compute_means(data_func, "indep_lag")
  compute_means(data_func, "indep")
  compute_demean(data_func, "indep_lag")
  compute_demean(data_func, "indep")

  # Define objective functions for Newton's method
  func_gamma <- function(rho) {
    sum((data_func$indep_demean - rho * data_func$indep_lag_demean)^2) / (num_individuals * (T0 - 2))
  }
  func_gamma_prime <- function(rho) {
    -2 * sum(data_func$indep_lag_demean * (data_func$indep_demean - rho * data_func$indep_lag_demean)) / (num_individuals * (T0 - 2))
  }

  func_g <- function(rho) {
    (1 - (1 - rho^(T0-1)) / ((T0-1) * (1 - rho))) / ((T0-1) * (1 - rho))
  }
  func_g_prime <- function(rho) {
    1 / ((T0 - 1) * (1 - rho)^2) + ((T0-1) * rho^(T0-2) * (1 - rho) - 2 * (1 - rho^(T0-1))) / ((T0-1)^2 * (1 - rho)^3)
  }

  func_f <- function(rho) {
    mean(data_func$indep_lag_demean * (data_func$indep_demean - rho * data_func$indep_lag_demean)) + func_gamma(rho) * func_g(rho)
  }
  func_f_prime <- function(rho) {
    mean(- data_func$indep_lag_demean^2) + func_gamma(rho) * func_g_prime(rho) + func_gamma_prime(rho) * func_g(rho)
  }

  Newton <- function(rho) {
    rho - func_f(rho) / func_f_prime(rho)
  }

  # Newton's method iteration
  rho_init <- 0.4
  gap <- 1
  epsilon <- 0.0001

  while (gap > epsilon) {
    rho_new <- Newton(rho_init)
    gap <- abs(rho_init - rho_new)
    rho_init <- rho_new
  }

  rho.hat <- rho_init
  s.u2 <- func_gamma(rho.hat)

  # Get initial coefficient estimate
  h = 0
  data_horizon <- prepare_data(data, h, lag = lag, dif = dif)
  data_demean <- demean_data(data_horizon, lag = lag, vars_num = vars_num)
  rm(data_horizon)
  gc()

  y.vec <- data_demean$dep_lead_h_demean
  x.vec <- data_demean$indep_demean
  all.fit <- summary(lm(y.vec~x.vec-1))
  b0.hat <- all.fit$coefficients[1]

  # Apply bias correction
  f.rho <- s.u2*( (1-rho.hat^(0:horizon)) - (0:horizon/(T0-0:horizon)))/(1-rho.hat)^2
  sx = sqrt(mean(x.vec^2))

  results_le[, beta := results_fe$beta + b0.hat * f.rho / ( (T0-0:horizon) * sx^2)]
  return(results_le)
}