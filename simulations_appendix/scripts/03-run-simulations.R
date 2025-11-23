# =============================================================================
# 03-run-simulations.R
# Purpose: Run Monte Carlo simulations for panel local projection estimators
# =============================================================================

# Enable global progress updates
handlers(global = TRUE)
handlers("progress")

# =============================================================================
# Single Iteration Function
# =============================================================================

# Run a single simulation iteration
run_single_iteration <- function(num_individuals = 30, num_periods = 60, ar_coefficient = 0.5,
                                 beta = -0.25, intercept = 0.2, fixed_effect_scale = 0.2,
                                 random_noise_sd = 1, common_shock_sd = 1, error_correlation = 0,
                                 lag_coefficient = 0.5, ar_lag_coefficient = -0.5,
                                 time_effects = TRUE, lag = FALSE, dif = FALSE,
                                 horizon = 10, significance_level = 0.05, cluster = FALSE,
                                 A = NULL, vars_num = 2,
                                 estimate_beta_se_list, methods, sigma = NULL,
                                 cluster_method = "driscoll_kraay") {
  # Generate synthetic data
  data <- generate_simulation_data(
    num_individuals = num_individuals,
    num_periods = num_periods,
    ar_coefficient = ar_coefficient,
    beta = beta,
    intercept = intercept,
    fixed_effect_scale = fixed_effect_scale,
    random_noise_sd = random_noise_sd,
    common_shock_sd = common_shock_sd,
    error_correlation = error_correlation,
    lag_coefficient = lag_coefficient,
    ar_lag_coefficient = ar_lag_coefficient,
    time_effects = time_effects,
    lag = lag,
    dif = dif,
    A = A,
    vars_num = vars_num,
    sigma = sigma
  )

  num_methods <- length(methods)

  if (vars_num == 2) {
    # Estimate using each method
    results <- vector("list", num_methods + 1)
    for (i in 1:num_methods) {
      results[[i]] <- estimate_beta_se_list[[methods[i]]](
        data, num_periods, num_individuals,
        horizon = horizon, time_effects = time_effects, lag = lag,
        cluster = cluster, dif = dif, vars_num = vars_num,
        cluster_method = cluster_method)
    }

    # Calculate true IRF values
    true_irfs <- beta * ar_coefficient^(0:horizon)

    # Adjust for lag case
    if (lag) {
      P <- matrix(
        c(
          lag_coefficient + beta * ar_lag_coefficient,
          ar_lag_coefficient,
          beta * ar_coefficient,
          ar_coefficient
        ),
        2, 2
      )
      P_h <- P
      true_irfs[1] <- beta
      true_irfs[2] <- beta * ar_coefficient
      for (h in 2:horizon) {
        P_h <- P_h %*% P
        true_irfs[h + 1] <- P_h[1, 2]
      }
    }

    # Add true IRF to results
    if (dif) {
      true_irfs <- beta * ar_coefficient * (1 - ar_coefficient^(1:horizon)) / (1 - ar_coefficient)
      results[[num_methods + 1]] <- data.table(beta = true_irfs, Horizon = 1:horizon, se = numeric(10), estimator = "true IRF")
    } else {
      results[[num_methods + 1]] <- data.table(beta = true_irfs, Horizon = 0:horizon, se = numeric(11), estimator = "true IRF")
    }

    # Combine results and compute metrics
    results_table <- rbindlist(results, use.names = TRUE)
    results_table[, true_beta := rep(true_irfs, times = num_methods + 1)]
    v <- qnorm(1 - significance_level / 2)
    results_table[, coverage := true_beta >= beta - v * se & true_beta <= beta + v * se]
    results_table[, rmse := (true_beta - beta)^2]
    results_table[, true_beta := NULL]
  } else {
    # 3-variable case: compute true IRFs
    true_irfsZ <- numeric(horizon)
    true_irfsY <- numeric(horizon)

    A_h <- A
    true_irfsZ[1] <- A[2, 1]
    true_irfsY[1] <- A[3, 1]
    for (h in 2:horizon) {
      A_h <- A_h %*% A
      true_irfsZ[h] <- A_h[2, 1]
      true_irfsY[h] <- A_h[3, 1]
    }

    # Estimate with Z as dependent variable
    setnames(data, c("x_var", "y_var", "z_var"), c("indep", "indep2", "dep"))
    resultsZ <- vector("list", num_methods + 1)
    for (i in 1:num_methods) {
      resultsZ[[i]] <- estimate_beta_se_list[[methods[i]]](
        data, num_periods, num_individuals,
        horizon = horizon, time_effects = time_effects, lag = lag,
        cluster = cluster, dif = dif, vars_num = vars_num,
        cluster_method = cluster_method)
    }
    resultsZ[[num_methods + 1]] <- data.table(beta = true_irfsZ, Horizon = 1:horizon, estimator = "true IRF")

    results_tableZ <- rbindlist(resultsZ, use.names = TRUE)
    results_tableZ[, true_beta := rep(true_irfsZ, times = num_methods + 1)]
    results_tableZ[, rmse := (true_beta - beta)^2]
    results_tableZ[, true_beta := NULL]
    results_tableZ[, dep := "Z"]

    # Estimate with Y as dependent variable
    setnames(data, c("indep", "indep2", "dep"), c("indep", "dep", "indep2"))
    resultsY <- vector("list", num_methods + 1)
    for (i in 1:num_methods) {
      resultsY[[i]] <- estimate_beta_se_list[[methods[i]]](
        data, num_periods, num_individuals,
        horizon = horizon, time_effects = time_effects, lag = lag,
        cluster = cluster, dif = dif, vars_num = vars_num,
        cluster_method = cluster_method)
    }
    resultsY[[num_methods + 1]] <- data.table(beta = true_irfsY, Horizon = 1:horizon, estimator = "true IRF")

    results_tableY <- rbindlist(resultsY, use.names = TRUE)
    results_tableY[, true_beta := rep(true_irfsY, times = num_methods + 1)]
    results_tableY[, rmse := (true_beta - beta)^2]
    results_tableY[, true_beta := NULL]
    results_tableY[, dep := "Y"]

    # Combine results and compute contrast
    results_table <- rbindlist(list(results_tableZ, results_tableY))
    results_table[, cumsum_beta := cumsum(beta), by = .(estimator, dep)]
    cumsum_Z = as.numeric(results_table[dep == "Z", cumsum_beta])
    cumsum_Y = as.numeric(results_table[dep == "Y", cumsum_beta])
    results_table[, cumsum_beta := NULL]
    results_table[, contrast_YZ := rep(cumsum_Y / cumsum_Z, 2)]
    results_table[, contrast_rmse := (contrast_YZ - contrast_YZ[estimator == "true IRF"])^2, by = .(dep)]
  }

  return(results_table)
}

# =============================================================================
# Parameter-Level Simulation Function
# =============================================================================

# Run simulations for a specific parameter combination
run_simulation_for_parameters <- function(num_individuals = 30, num_periods = 60, ar_coefficient = 0.5,
                                          beta = -0.25, intercept = 0.2, fixed_effect_scale = 0.2,
                                          random_noise_sd = 1, common_shock_sd = 1, error_correlation = 0,
                                          lag_coefficient = 0.5, ar_lag_coefficient = -0.5,
                                          time_effects = TRUE, lag = FALSE, dif = FALSE,
                                          horizon = 10, significance_level = 0.05, cluster = FALSE,
                                          estimate_beta_se_list, methods, A = NULL, vars_num = 2,
                                          num_iterations = 100, runplace = "local", parallel = FALSE, sigma = NULL,
                                          cluster_method = "driscoll_kraay") {

  if (parallel) {
    # Parallel execution
    if (runplace == "cloud") {
      results <- with_progress({
        p <- progressr::progressor(along = seq_len(num_iterations))
        future_map(seq_len(num_iterations), function(iter) {
          p(sprintf("Iteration %d/%d", iter, num_iterations))
          source("scripts/01-generate-simulation-data.R")
          source("scripts/02-estimate.R")
          run_single_iteration(
            num_individuals = num_individuals, num_periods = num_periods,
            ar_coefficient = ar_coefficient, beta = beta, intercept = intercept,
            fixed_effect_scale = fixed_effect_scale, random_noise_sd = random_noise_sd,
            common_shock_sd = common_shock_sd, error_correlation = error_correlation,
            lag_coefficient = lag_coefficient, ar_lag_coefficient = ar_lag_coefficient,
            time_effects = time_effects, lag = lag, dif = dif,
            horizon = horizon, significance_level = significance_level, cluster = cluster,
            estimate_beta_se_list = estimate_beta_se_list, methods = methods,
            A = A, vars_num = vars_num, sigma = sigma,
            cluster_method = cluster_method
          )
        }, .options = furrr_options(seed = TRUE))
      })
    } else {
      handlers("rstudio", "progress")
      results <- with_progress({
        p <- progressr::progressor(along = seq_len(num_iterations))
        future_map(seq_len(num_iterations), function(iter) {
          p(sprintf("Iteration %d/%d", iter, num_iterations))
          source("scripts/01-generate-simulation-data.R")
          source("scripts/02-estimate.R")
          run_single_iteration(
            num_individuals = num_individuals, num_periods = num_periods,
            ar_coefficient = ar_coefficient, beta = beta, intercept = intercept,
            fixed_effect_scale = fixed_effect_scale, random_noise_sd = random_noise_sd,
            common_shock_sd = common_shock_sd, error_correlation = error_correlation,
            lag_coefficient = lag_coefficient, ar_lag_coefficient = ar_lag_coefficient,
            time_effects = time_effects, lag = lag, dif = dif,
            horizon = horizon, significance_level = significance_level, cluster = cluster,
            estimate_beta_se_list = estimate_beta_se_list, methods = methods,
            A = A, vars_num = vars_num, sigma = sigma,
            cluster_method = cluster_method
          )
        }, .options = furrr_options(seed = TRUE))
      })
    }
  } else {
    # Sequential execution
    if (runplace == "cloud") {
      pboptions(type = "timer")
      pb <- txtProgressBar(min = 0, max = num_iterations, style = 3)
      results <- pblapply(seq_len(num_iterations), function(iter) {
        setTxtProgressBar(pb, iter)
        flush.console()
        run_single_iteration(
          num_individuals = num_individuals, num_periods = num_periods,
          ar_coefficient = ar_coefficient, beta = beta, intercept = intercept,
          fixed_effect_scale = fixed_effect_scale, random_noise_sd = random_noise_sd,
          common_shock_sd = common_shock_sd, error_correlation = error_correlation,
          lag_coefficient = lag_coefficient, ar_lag_coefficient = ar_lag_coefficient,
          time_effects = time_effects, lag = lag, dif = dif,
          horizon = horizon, significance_level = significance_level, cluster = cluster,
          estimate_beta_se_list = estimate_beta_se_list, methods = methods,
          A = A, vars_num = vars_num, sigma = sigma,
          cluster_method = cluster_method
        )
      })
    } else {
      results <- pblapply(seq_len(num_iterations), function(iter) {
        run_single_iteration(
          num_individuals = num_individuals, num_periods = num_periods,
          ar_coefficient = ar_coefficient, beta = beta, intercept = intercept,
          fixed_effect_scale = fixed_effect_scale, random_noise_sd = random_noise_sd,
          common_shock_sd = common_shock_sd, error_correlation = error_correlation,
          lag_coefficient = lag_coefficient, ar_lag_coefficient = ar_lag_coefficient,
          time_effects = time_effects, lag = lag, dif = dif,
          horizon = horizon, significance_level = significance_level, cluster = cluster,
          estimate_beta_se_list = estimate_beta_se_list, methods = methods,
          A = A, vars_num = vars_num, sigma = sigma,
          cluster_method = cluster_method
        )
      })
    }
  }

  # Aggregate results
  conbined_results <- rbindlist(results, use.names = TRUE)
  if (vars_num == 2) {
    mean_results <- conbined_results[, lapply(.SD, mean), by = .(estimator, Horizon)]
  } else {
    mean_results <- conbined_results[, lapply(.SD, mean), by = .(estimator, Horizon, dep)]
    mean_results[, contrast_rmse := sqrt(contrast_rmse)]
  }
  mean_results[, rmse := sqrt(rmse)]

  return(mean_results)
}

# =============================================================================
# Main Simulation Function
# =============================================================================

# Run simulations across multiple parameter combinations
run_simulations <- function(beta = -0.25, intercept = 0.2, fixed_effect_scale = 0.2,
                            random_noise_sd = 1, common_shock_sd = 1, error_correlation = 0,
                            lag_coefficient = 0.5, ar_lag_coefficient = -0.5,
                            time_effects = TRUE, lag = FALSE, dif = FALSE,
                            horizon = 10, significance_level = 0.05, cluster = FALSE,
                            estimate_beta_se_list, methods, num_iterations = 100,
                            A = NULL, vars_num = 2,
                            ar_coefficients, parameter_combinations, runplace = "local", parallel = FALSE,
                            cluster_method = "driscoll_kraay") {

  if (vars_num == 2) {
    # 2-variable case
    num_situations <- length(ar_coefficients) * length(parameter_combinations)
    results <- vector("list", num_situations)
    i <- 1

    for (params in parameter_combinations) {
      num_individuals <- params[1]
      num_periods <- params[2]

      # Generate covariance matrix for correlated errors
      if (error_correlation > 0) {
        lambda <- runif(num_individuals, min = 0, max = 1 / sqrt(2))
        sigma <- lambda %*% t(lambda)
        diag(sigma) <- 1
      } else {
        sigma = 0
      }

      for (ar_coeff in ar_coefficients) {
        message("Simulating for AR: ", ar_coeff, " N: ", num_individuals, " T: ", num_periods)

        results[[i]] <- run_simulation_for_parameters(
          num_individuals = num_individuals, num_periods = num_periods,
          ar_coefficient = ar_coeff, beta = beta, intercept = intercept,
          fixed_effect_scale = fixed_effect_scale, random_noise_sd = random_noise_sd,
          common_shock_sd = common_shock_sd, error_correlation = error_correlation,
          lag_coefficient = lag_coefficient, ar_lag_coefficient = ar_lag_coefficient,
          time_effects = time_effects, lag = lag, dif = dif,
          horizon = horizon, significance_level = significance_level, cluster = cluster,
          estimate_beta_se_list = estimate_beta_se_list, methods = methods,
          num_iterations = num_iterations, runplace = runplace, parallel = parallel,
          A = A, vars_num = vars_num, sigma = sigma,
          cluster_method = cluster_method
        )

        results[[i]][, `:=`(N = num_individuals, T = num_periods, AR_Coefficient = ar_coeff)]
        i <- i + 1
      }
    }

  } else {
    # 3-variable case
    num_situations <- length(parameter_combinations)
    results <- vector("list", num_situations)
    i <- 1

    for (params in parameter_combinations) {
      num_individuals <- params[1]
      num_periods <- params[2]
      message("Simulating for N: ", num_individuals, " T: ", num_periods)

      results[[i]] <- run_simulation_for_parameters(
        num_individuals = num_individuals, num_periods = num_periods,
        ar_coefficient = 0, beta = beta, intercept = intercept,
        fixed_effect_scale = fixed_effect_scale, random_noise_sd = random_noise_sd,
        common_shock_sd = common_shock_sd, error_correlation = error_correlation,
        lag_coefficient = lag_coefficient, ar_lag_coefficient = ar_lag_coefficient,
        time_effects = time_effects, lag = lag, dif = dif,
        horizon = horizon, significance_level = significance_level, cluster = cluster,
        estimate_beta_se_list = estimate_beta_se_list, methods = methods,
        num_iterations = num_iterations, runplace = runplace, parallel = parallel,
        A = A, vars_num = vars_num,
        cluster_method = cluster_method
      )

      results[[i]][, `:=`(N = num_individuals, T = num_periods)]
      i <- i + 1
    }
  }

  return(rbindlist(results, use.names = TRUE))
}
