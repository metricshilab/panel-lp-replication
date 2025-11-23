# =============================================================================
# 02-01-data-preparation.R
# Purpose: Prepare and transform data for panel local projection estimation
# =============================================================================

# =============================================================================
# Data Preparation Functions
# =============================================================================

# Prepare data for a given horizon
prepare_data <- function(data, horizon, lag = FALSE, dif = FALSE) {
  data <- copy(data)

  # Create lead response variable
  data[, dep_lead_h := shift(dep, horizon, type = "lead"), by = id]
  data <- na.omit(data, cols = "dep_lead_h")

  # Apply differencing if specified
  if (dif) {
    data[,dep_lead_h := dep_lead_h - dep, by = id]
  }

  # Add lagged dependent variable if specified
  if (lag) {
    if (horizon >= 1) {
      data[, dep_lag := dep, by =id]
    } else {
      data[, dep_lag := shift(dep, 1, type = "lag"), by = id]
    }
    data <- na.omit(data, cols = "dep_lag")
  }

  data[, dep := NULL]
  return(data)
}

# =============================================================================
# Mean Computation Functions
# =============================================================================

# Compute means at individual and time levels
compute_means <- function(data, var, midpoint = NULL) {
  # Individual and time-level means
  data[, paste0(var, "_time") := mean(get(var), na.rm = TRUE), by = time]
  data[, paste0(var, "_individual") := mean(get(var), na.rm = TRUE), by = id]
  data[, paste0(var, "_mean") := mean(get(var), na.rm = TRUE)]

  # Midpoint-based means for split-panel jackknife
  if (!is.null(midpoint)) {
    data[, paste0(var, "_individual_ab") :=
           ifelse(time <= midpoint,
                  mean(get(var)[time <= midpoint], na.rm = TRUE),
                  mean(get(var)[time > midpoint], na.rm = TRUE)),
         by = id]

    data[, paste0(var, "_mean_ab") :=
           ifelse(time <= midpoint,
                  mean(get(var)[time <= midpoint], na.rm = TRUE),
                  mean(get(var)[time > midpoint], na.rm = TRUE))]
  }
}

# =============================================================================
# Demeaning Functions
# =============================================================================

# Compute demeaned values
compute_demean <- function(data, var, midpoint = NULL) {
  # Standard demeaning
  data[, paste0(var, "_demean") := get(var) - get(paste0(var, "_individual")) -
         get(paste0(var, "_time")) + get(paste0(var, "_mean"))]

  # Midpoint-adjusted demeaning for split-panel jackknife
  if (!is.null(midpoint)) {
    data[, paste0(var, "_demean_ab") := get(var) - get(paste0(var, "_individual_ab")) -
           get(paste0(var, "_time")) + get(paste0(var, "_mean_ab"))]
  }
}

# Demean data by group with time effects
demean_data <- function(data, midpoint = NULL, lag = FALSE, vars_num = 2) {
  data <- copy(data)

  # Compute means and demeaned values for main variables
  compute_means(data, "dep_lead_h", midpoint)
  compute_means(data, "indep", midpoint)
  compute_demean(data, "dep_lead_h", midpoint)
  compute_demean(data, "indep", midpoint)

  # Handle lagged dependent variable
  if (lag) {
    compute_means(data, "dep_lag", midpoint)
    compute_demean(data, "dep_lag", midpoint)
  }

  # Handle second independent variable for 3-variable case
  if (vars_num == 3) {
    compute_means(data, "indep2", midpoint)
    compute_demean(data, "indep2", midpoint)
  }

  # Select relevant columns
  selected_cols <- grep("demean$|demean_ab$|^id$|^time$|^indep$|^indep_individual_ab$|^indep2_individual_ab$|^dep_lag_individual_ab$|^dep_lag$|^indep2$", names(data), value = TRUE)
  data <- data[, ..selected_cols]

  return(data)
}

# =============================================================================
# DGMM Differencing Functions
# =============================================================================

# Compute differenced values for DGMM
compute_diff_dgmm <- function(data, var) {
  data[, paste0(var, "_demean") := mean(get(var), na.rm = TRUE) - get(var), by = time]
  data[, paste0(var, "_diff") := c(NA, diff(get(paste0(var, "_demean")))), by = id]
}

# Difference data with time effects for DGMM estimation
difference_data <- function(data) {
  data <- copy(data)

  # Demean and difference variables
  compute_diff_dgmm(data, "dep_lag")
  compute_diff_dgmm(data, "indep")
  compute_diff_dgmm(data, "dep_lead_h")

  # Create instrumental variables
  data[, iv_dep := shift(dep_lag_demean, 1, type = "lag"), by = id]
  data[, iv_indep := shift(indep_demean, 1, type = "lag"), by = id]

  data <- na.omit(data)
  return(data)
}
