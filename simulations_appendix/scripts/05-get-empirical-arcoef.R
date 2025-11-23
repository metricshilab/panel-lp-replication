# =============================================================================
# 05-get-empirical-arcoef.R
# Purpose: Estimate empirical AR(1) coefficient from currency crisis data
# =============================================================================

# Load required libraries
library(tidyverse)
library(readxl)
library(pLP)

# =============================================================================
# Main Function
# =============================================================================

# Estimate AR(1) coefficient using Split-Panel Jackknife method
get_empirical_arcoef <- function(){

  # Load and prepare data
  data <- read_excel("data/currcrisis/20050666_data.xls",
                     sheet = "currcrisis") |>
    filter(obs >= 1965 & obs <= 2000) |>
    select(where(~ !any(is.na(.)))) |>
    pivot_longer(-obs, names_to = "country", values_to = "y") |>
    group_by(country) |>
    mutate(y_lag = lag(y)) |>
    ungroup() |>
    drop_na()

  # Rename time variable
  data <- data |>
    rename(time = obs)

  # Generate numeric ID for each country
  data$id <- as.integer(factor(data$country, levels = unique(data$country)))
  data <- data |> select(id, everything())

  # Ensure correct data types
  data$id <- as.integer(data$id)
  data$time <- as.integer(data$time)

  # Sort and select required columns
  data <- data[order(data$id, data$time), ]
  data <- data[, c("id", "time", "y", "y_lag")]
  df <- as.data.frame(data)

  # Estimate AR coefficient using panel local projection
  results <- panelLP(
    df, Y.name = c("y"), X.name = c("y_lag"),
    method = "SPJ", te = FALSE, H = 0)

  return(as.numeric(results$IRF))
}
