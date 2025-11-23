# =============================================================================
# 01-generate-simulation-data.R
# Purpose: Generate synthetic panel data for Monte Carlo simulations
# =============================================================================

# =============================================================================
# Helper Functions
# =============================================================================

# Generate individual-level fixed effects
generate_fixed_effects <- function(indep_means, fixed_effect_scale, common_shock_sd, num_periods) {
  common_shock <- rnorm(1, 0, common_shock_sd)
  return(fixed_effect_scale * sqrt(num_periods) * indep_means + common_shock)
}

# Generate independent variable as AR(1) process
generate_indep <- function(num_periods, ar_coefficient, intercept, random_noise) {
  x_vals <- numeric(num_periods)
  x_vals[1] <- intercept + ar_coefficient * rnorm(1) + random_noise[1]
  for (t in 2:num_periods) {
    x_vals[t] <- intercept + ar_coefficient * x_vals[t - 1] + random_noise[t]
  }
  return(x_vals)
}

# Generate dependent variable
generate_dep <- function(num_periods, fixed_effect, indep, beta, error_term, time_effects) {
  y_vals <- numeric(num_periods)
  for (t in 1:num_periods) {
    y_vals[t] <- fixed_effect[t] + beta * indep[t] + error_term[t] +
      ifelse(time_effects, (0.025 * t + 0.001 * t^2), 0)
  }
  return(y_vals)
}

# Generate correlated error terms with cross-sectional dependence
generate_error_terms <- function(num_individuals, correlation, sigma) {
  lambda <- runif(num_individuals, min = 0, max = 1 / sqrt(2))
  sigma <- lambda %*% t(lambda)
  diag(sigma) <- 1
  return(as.vector(mvtnorm::rmvnorm(n = 1, mean = rep(0, num_individuals), sigma = sigma)))
}

# =============================================================================
# Main Data Generation Functions
# =============================================================================

# Dispatcher function to generate simulation data
generate_simulation_data <- function(num_individuals = 30, num_periods = 60, ar_coefficient = 0.5,
                                           beta = -0.25, intercept = 0.2, fixed_effect_scale = 0.2,
                                           random_noise_sd = 1, common_shock_sd = 1, error_correlation = 0,
                                           lag_coefficient = 0.5, ar_lag_coefficient = -0.5,
                                           time_effects = TRUE, lag = FALSE, dif = FALSE,
                                           A = NULL, vars_num = 2, sigma = NULL) {
  if (vars_num == 2) {
    return(generate_simulation_data_vars2(num_individuals, num_periods, ar_coefficient, beta, intercept,
                                          fixed_effect_scale, random_noise_sd, common_shock_sd, error_correlation,
                                          lag_coefficient, ar_lag_coefficient, time_effects, lag, dif, sigma))
  } else if (vars_num == 3) {
    return(generate_simulation_data_vars3(num_individuals, num_periods, A))
  }
}

# Generate simulation data for 2-variable case
generate_simulation_data_vars2 <- function(num_individuals = 30, num_periods = 60, ar_coefficient = 0.5,
                                     beta = -0.25, intercept = 0.2, fixed_effect_scale = 0.2,
                                     random_noise_sd = 1, common_shock_sd = 1, error_correlation = 0,
                                     lag_coefficient = 0.5, ar_lag_coefficient = -0.5,
                                     time_effects = TRUE, lag = FALSE, dif = FALSE, sigma = NULL) {

  # Create base data.table
  dt <- CJ(id = 1:num_individuals, time = 1:num_periods)

  # Generate random noise
  dt[, random_noise := rnorm(.N, 0, random_noise_sd)]

  # Generate error terms
  if (error_correlation == 0) {
    dt[, error_term := rnorm(.N)]
  } else {
    error_terms <- replicate(num_periods, generate_error_terms(num_individuals, error_correlation, sigma), simplify = TRUE)
    error_terms_vec <- as.vector(t(error_terms))
    dt[, error_term := error_terms_vec]
  }
  
  # Generate variables based on lag specification
  if (lag) {
    # Initialize variables
    dt[, indep := 0.0]
    dt[time == 1, indep := intercept + ar_coefficient * rnorm(num_individuals, 0, 1) + random_noise]

    dt[, dep := 0.0]
    dt[, fixed_effect := generate_fixed_effects(0, fixed_effect_scale, common_shock_sd, num_periods), by = id]
    dt[time == 1, dep := fixed_effect + beta * indep + error_term + ifelse(time_effects, (0.025 + 0.001), 0), by = id]

    # Iterate through time periods
    for (t in 2:num_periods) {
      dt[time == t, c("indep", "dep") := {
        prev <- dt[time == t - 1, .(indep, dep), on = "id"]
        prev <- t(as.matrix(prev))
        new_values <- matrix(0, nrow = 2, ncol = .N)
        new_values[1,] <- intercept + ar_lag_coefficient * prev[2,] + ar_coefficient * prev[1,] + random_noise
        new_values[2,] <- fixed_effect + beta * new_values[1,] + error_term + lag_coefficient * prev[2,] +
          ifelse(time_effects, (0.025 * t + 0.001 * t^2), 0)

        as.data.table(t(new_values))
      }]
    }
  } else {
    # Generate independent variable
    dt[, indep := generate_indep(num_periods, ar_coefficient, intercept, random_noise), by = id]

    # Compute means and fixed effects
    dt[, indep_mean := mean(indep), by = id]
    dt[, fixed_effect := generate_fixed_effects(indep_mean, fixed_effect_scale, common_shock_sd, num_periods), by = id]

    # Generate dependent variable
    if (dif) {
      dt[, dif_dep := generate_dep(num_periods, fixed_effect, indep, beta, error_term, time_effects), by = id]
      dt[, dep := cumsum(dif_dep), by = id]
    } else {
      dt[, dep := generate_dep(num_periods, fixed_effect, indep, beta, error_term, time_effects), by = id]
    }
  }

  # Return selected columns
  dt <- dt[, .(id, time, dep, indep)]
  return(dt)
}

# Generate simulation data for 3-variable VAR case
generate_simulation_data_vars3 <- function(num_individuals, num_periods, A) {
  # Create data.table
  dt <- data.table(expand.grid(id = 1:num_individuals, time = 0:num_periods))

  # Generate fixed effects
  fixed_effects <- data.table(
    id = 1:num_individuals,
    alpha_X = rnorm(num_individuals, 0, 1),
    alpha_Z = rnorm(num_individuals, 0, 1),
    alpha_Y = rnorm(num_individuals, 0, 1)
  )
  dt <- merge(dt, fixed_effects, by = "id", all.x = TRUE)

  # Initialize variables
  dt[, `:=`(x_var = 0, z_var = 0, y_var = 0)]

  # Set initial values at t = 0
  dt[time == 0, `:=`(
    x_var = alpha_X + rnorm(.N, 0, 1),
    z_var = alpha_Z + rnorm(.N, 0, 1),
    y_var = alpha_Y + rnorm(.N, 0, 1)
  )]

  # Simulate VAR process
  for (t in 1:num_periods) {
    dt[time == t, c("x_var", "z_var", "y_var") := {
      prev <- dt[time == t - 1, .(x_var, z_var, y_var), on = "id"]
      epsilon <- matrix(rnorm(3 * .N), ncol = 3)
      new_values <- A %*% t(as.matrix(prev)) + t(epsilon)

      # Add fixed effects
      new_values[1, ] <- new_values[1, ] + alpha_X
      new_values[2, ] <- new_values[2, ] + alpha_Z
      new_values[3, ] <- new_values[3, ] + alpha_Y

      as.data.table(t(new_values))
    }]
  }

  dt[time > 0, .(id, time, x_var, z_var, y_var)]
}