runplace = "local" # pick from c("local", "cloud")
parallel = TRUE
cluster_num = 64
# to monitor the progress in the cloud, type this in the terminal
# export R_PROGRESSR_ENABLE=TRUE

start_time <- Sys.time()
message("Execution started at:", start_time)

packages = c("mvtnorm", "data.table", "AER", "future", "furrr", "progressr", "pbapply")
# Load all the required packages
invisible(lapply(packages, require, character.only = TRUE))
# Functions -----------------------------------------------------------------

source("scripts/01-generate-simulation-data.R")
source("scripts/02-estimate.R")
source("scripts/05-get-empirical-arcoef.R")
if (runplace == "local") {
  # load packages for plot
  packages_plot = c("ggplot2", "tidyverse", "extrafont")
  invisible(lapply(packages_plot, require, character.only = TRUE))
  source("scripts/04-plot.R")
}
source("scripts/03-run-simulations.R")

# global debugging
options(error = recover)
# Parameters --------------------------------------------------------------
source("config.R")

# Set cases and methods
cases <- config$cases
methods <- config$methods
twoway_cluster <- config$twoway_cluster
cluster_method <- config$cluster_method

# DGP parameters
beta <- config$beta
ar_coefficients <- config$ar_coefficients
intercept <- 0.2
fixed_effect_scale <- 0.2
random_noise_sd <- 1
common_shock_sd <- 1
time_effects <- list(
  correlation = TRUE, 
  lag = TRUE, 
  highrho_simple = FALSE, 
  dif = FALSE, 
  simple = FALSE
)
dif <- list(
  correlation = FALSE, 
  lag = FALSE, 
  highrho_simple = FALSE, 
  dif = TRUE,
  simple = FALSE
)

# lag parameters
lag <- list(
  correlation = FALSE, 
  lag = TRUE, 
  highrho_simple = FALSE, 
  dif = FALSE, 
  simple = FALSE
)
lag_coefficient <- 0.5
ar_lag_coefficient <- -0.5

# the error_correlation of y individuals. represent cross-sectional correlation.
error_correlation <- list(
  correlation = 0.5, 
  lag = 0, 
  highrho_simple = 0, 
  dif = 0, 
  simple = 0
)

# vars3 case
A <- config$A
vars_num <- list(
  correlation = 2, 
  lag = 2, 
  highrho_simple = 2, 
  dif = 2, 
  simple = 2
)


# Estimation methods
estimate_beta_se_list <- list(
  spj = estimate_spj,
  fe = estimate_fe,
  dgmm = estimate_dgmm,
  db = estimate_db,
  le = estimate_le
)


# Simulation parameters
significance_level <- 0.05
horizon <- 10
num_iterations <- config$num_iterations
parameter_combinations <- config$parameter_combinations

# Simulation --------------------------------------------------------------

# Set seed
set.seed(2023)

# Set up the parallel processing plan to utilize multiple CPU cores
if (runplace == "cloud") {
  plan(multisession, workers = cluster_num) # in SCRP
} else {
  plan(multisession, workers = parallel::detectCores() - 1) # in local computer
}

# plan(sequential)

for (case in cases) {
  for (cluster in twoway_cluster) {
    message("Run ", num_iterations, " simulations for case: ", case)
    simulation_results <- run_simulations(
      beta = beta[[case]],
      intercept = intercept,
      fixed_effect_scale = fixed_effect_scale,
      random_noise_sd = random_noise_sd,
      common_shock_sd = common_shock_sd,
      error_correlation = error_correlation[[case]],
      lag_coefficient = lag_coefficient,
      ar_lag_coefficient = ar_lag_coefficient,
      time_effects = time_effects[[case]],
      lag = lag[[case]],
      dif = dif[[case]],
      horizon = horizon,
      significance_level = significance_level,
      cluster = cluster,
      estimate_beta_se_list = estimate_beta_se_list,
      methods = methods,
      num_iterations = num_iterations,
      ar_coefficients = ar_coefficients[[case]],
      parameter_combinations = parameter_combinations[[case]],
      runplace = runplace,
      parallel = parallel,
      A = A[[case]],
      vars_num = vars_num[[case]],
      cluster_method = cluster_method
    )
    
    # Save simulation results to a file
    write_csv(
      simulation_results, 
      paste0("results/tables/results_", case, "_", 
             ifelse(cluster, "twoway_se", "simple_se"), ".csv")
    )
  }
}



# Plotting ----------------------------------------------------------------

# Parameters of the plotting
limits_vs <- list(
  correlation = c(-0.7, 0.1), 
  lag = c(-0.3, 0.1), 
  highrho_simple = c(-0.7, 0.1),
  dif = c(-0.2, 0.2),
  simple = c(-0.7, 0.1)
)
breaks_vs <- list(
  correlation = c(-0.6, -0.4, -0.2, 0.0), 
  lag = c(-0.3, -0.2, -0.1, 0.0), 
  highrho_simple = c(-0.6, -0.4, -0.2, 0.0),
  dif = c(-0.2, -0.1, 0.0, 0.1),
  simple = c(-0.6, -0.4, -0.2, 0.0)
)

# Prepare results for plotting
for (case in cases) {
  for (cluster in twoway_cluster){
    results <- read_csv(paste0("results/tables/results_", case, "_",
                               ifelse(cluster, "twoway_se", "simple_se"), ".csv"))
  }
  # Add confidence interval length
  results <- results |>
    mutate(
      confidence_interval_length = 2 * qnorm(1 - significance_level / 2) * se
    )
  
  plot_simulation_results(
    results,
    limits_vs = limits_vs[[case]],
    breaks_vs = breaks_vs[[case]],
    case = case
  )
}

end_time <- Sys.time()
message("Execution ended at:", end_time)
# Calculate and display the time difference
time_used <- end_time - start_time
message("Time used:", time_used)
