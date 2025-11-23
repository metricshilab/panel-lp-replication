# =============================================================================
# 02-03-error-calculation.R
# Purpose: Calculate standard errors for panel local projection estimators
# =============================================================================

# =============================================================================
# Data Preparation for Standard Errors
# =============================================================================

# Prepare data for SPJ standard error calculation
longer_data_se_spj <- function(data, midpoint, beta, lag = FALSE) {
  data <- copy(data)
  data[, d := 2 * indep_demean - indep_demean_ab]

  if (lag) {
    data[, d_lag := 2 * dep_lag_demean - dep_lag_demean_ab]
    data[, residual := dep_lead_h_demean - beta[1] * indep_demean - beta[2] * dep_lag_demean]
    data[, d_lag_e := d_lag * residual]
  } else {
    data[, residual := dep_lead_h_demean - beta * indep_demean]
  }

  data[, d_e := d * residual]
  selected_cols <- grep("^d_e$|^id$|^time$|^indep_demean$|^dep_lag_demean$|^d_lag_e$", names(data), value = TRUE)
  data <- data[, ..selected_cols]

  return(data)
}

# Prepare data for FE standard error calculation
longer_data_se_fe <- function(data, beta, lag = FALSE) {
  data <- copy(data)

  if (lag) {
    data[, residual := dep_lead_h_demean - beta[1] * indep_demean - beta[2] * dep_lag_demean]
    data[, `:=`(
      d_e = indep_demean * residual,
      d_lag_e = dep_lag_demean * residual
    )]
  } else {
    data[, residual := dep_lead_h_demean - beta * indep_demean]
    data[, d_e := indep_demean * residual]
  }

  selected_cols <- grep("^d_e$|^id$|^time$|^indep_demean$|^dep_lag_demean$|^d_lag_e$", names(data), value = TRUE)
  data <- data[, ..selected_cols]

  return(data)
}

# Prepare data for DGMM standard error calculation
longer_data_se_dgmm <- function(data, residuals) {
  data <- copy(data)
  data[, residual := unlist(residuals)]
  data[, `:=`(d_e = iv_indep * residual, d_lag_e = iv_dep * residual)]
  return(data)
}

# =============================================================================
# Q Matrix Calculation
# =============================================================================

# Calculate Q matrix for variance estimation
calculate_Q <- function(data, N, lag) {
  T_h <- uniqueN(data$time)

  if (lag) {
    X <- as.matrix(data[, .(indep_demean, dep_lag_demean)])
  } else {
    X <- as.matrix(data[, .(indep_demean)])
  }

  Q <- t(X) %*% X / (N * T_h)
  return(Q)
}

# =============================================================================
# SPJ Standard Error Calculation
# =============================================================================

# Calculate standard errors for SPJ estimates
calculate_se_spj <- function(data, N, midpoint, beta, lag = FALSE, cluster = FALSE,
                             cluster_method = "driscoll_kraay") {
  T_h <- uniqueN(data$time)

  # Compute Q matrix from both periods
  data_a <- data[time <= midpoint]
  data_b <- data[time > midpoint]
  Q_a <- calculate_Q(data_a, N, lag)
  Q_b <- calculate_Q(data_b, N, lag)
  Q <- 0.5 * (Q_a + Q_b)

  if (lag) {
    # Compute R_N for multivariate case
    N_long <- data[, c(rbind(d_e, d_lag_e)), by = id]
    N_matrix <- matrix(N_long$V1, ncol = N)
    N_result <- N_matrix %*% t(N_matrix)

    # Extract 2x2 submatrices
    N_indices <- seq(1, 2 * T_h, by = 2)
    N_index_grid <- expand.grid(row_start = N_indices, col_start = N_indices)
    N_submatrix <- apply(N_index_grid, 1, function(idx) {
      row <- idx[1]
      col <- idx[2]
      N_result[row:(row + 1), col:(col + 1)]
    })
    N_sum <- apply(N_submatrix, 1, sum)
    R_N <- matrix(N_sum, nrow = 2) / (N * T_h)

    # Compute variance
    if (cluster) {
      # Two-way clustered: compute R_T
      T_long <- data[, c(rbind(d_e, d_lag_e)), by = time]
      T_matrix <- matrix(T_long$V1, ncol = T_h)
      T_result <- T_matrix %*% t(T_matrix)

      T_indices <- seq(1, 2 * N, by = 2)
      T_index_grid <- expand.grid(row_start = T_indices, col_start = T_indices)
      T_submatrix <- apply(T_index_grid, 1, function(idx) {
        row <- idx[1]
        col <- idx[2]
        T_result[row:(row + 1), col:(col + 1)]
      })
      T_sum <- apply(T_submatrix, 1, sum)
      R_T <- matrix(T_sum, nrow = 2) / (N * T_h)

      # Compute R_NT
      NT_de <- as.matrix(data[, .(d_e, d_lag_e)])
      R_NT <- t(NT_de) %*% NT_de / (N * T_h)

      V_TW <- solve(Q) %*% (R_N + R_T - R_NT) %*% solve(Q)
    } else {
      V_TW <- solve(Q) %*% (R_N) %*% solve(Q)
    }

    ifelse(V_TW[1, 1] < 0, 0, sqrt(V_TW[1, 1] / (N * T_h)))

  } else {
    # Compute R_N for scalar case
    N_long <- data[, c(rbind(d_e)), by = id]
    N_matrix <- matrix(N_long$V1, ncol = N)
    N_result <- N_matrix %*% t(N_matrix)
    N_sum <- sum(N_result)
    R_N <- N_sum / (N * T_h)

    # Compute variance
    if (cluster) {
      if (cluster_method == "driscoll_kraay") {
        # Driscoll and Kraay (1998) standard errors
        g <- data[, sum(d_e), by = time]$V1
        m_T <- floor(T_h^(1/4))

        # Helper function for autocovariance calculation
        calculate_Delta_j <- function(g, T_h, j) {
          valid_indices <- (j + 1):T_h
          g_t <- g[valid_indices]
          g_t_minus_j <- g[valid_indices - j]
          Delta_j <- sum(g_t * g_t_minus_j) / T_h
          return(Delta_j)
        }

        # HAC variance estimation
        Delta_0 <- calculate_Delta_j(g, T_h, j = 0)
        R_NT_DK <- Delta_0

        if (m_T >= 1) {
          w <- 1 - 1:m_T / (m_T + 1)
          for (j in 1:m_T) {
            Delta_j <- calculate_Delta_j(g, T_h, j)
            R_NT_DK <- R_NT_DK + 2 * w[j] * Delta_j
          }
        }

        R_NT_DK <- R_NT_DK / (N * T_h)
        V_TW <- solve(Q) %*% (R_N + R_NT_DK) %*% solve(Q)

      } else {
        # Two-way clustered standard errors
        # Compute R_T
        T_long <- data[, c(rbind(d_e)), by = time]
        T_matrix <- matrix(T_long$V1, ncol = T_h)
        T_result <- T_matrix %*% t(T_matrix)
        T_sum <- sum(T_result)
        R_T <- T_sum / (N * T_h)

        # Compute R_NT
        NT_de <- as.matrix(data[, .(d_e)])
        R_NT <- t(NT_de) %*% NT_de / (N * T_h)

        V_TW <- solve(Q) %*% (R_N + R_T - R_NT) %*% solve(Q)
      }

    } else {
      V_TW <- solve(Q) %*% (R_N) %*% solve(Q)
    }

    return(as.numeric(ifelse(V_TW < 0, 0, sqrt(V_TW / (N * T_h)))))
  }
}

# =============================================================================
# FE Standard Error Calculation
# =============================================================================

# Calculate standard errors for FE estimates
calculate_se_fe <- function(data, N, lag = FALSE, cluster = FALSE,
                            cluster_method = "driscoll_kraay") {
  T_h <- uniqueN(data$time)
  Q <- calculate_Q(data, N, lag)

  if (lag) {
    # Compute R_N for multivariate case
    N_long <- data[, c(rbind(d_e, d_lag_e)), by = id]
    N_matrix <- matrix(N_long$V1, ncol = N)
    N_result <- N_matrix %*% t(N_matrix)

    # Extract 2x2 submatrices
    N_indices <- seq(1, 2 * T_h, by = 2)
    N_index_grid <- expand.grid(row_start = N_indices, col_start = N_indices)
    N_submatrix <- apply(N_index_grid, 1, function(idx) {
      row <- idx[1]
      col <- idx[2]
      N_result[row:(row + 1), col:(col + 1)]
    })
    N_sum <- apply(N_submatrix, 1, sum)
    R_N <- matrix(N_sum, nrow = 2) / (N * T_h)

    # Compute variance
    if (cluster) {
      # Two-way clustered: compute R_T
      T_long <- data[, c(rbind(d_e, d_lag_e)), by = time]
      T_matrix <- matrix(T_long$V1, ncol = T_h)
      T_result <- T_matrix %*% t(T_matrix)

      T_indices <- seq(1, 2 * N, by = 2)
      T_index_grid <- expand.grid(row_start = T_indices, col_start = T_indices)
      T_submatrix <- apply(T_index_grid, 1, function(idx) {
        row <- idx[1]
        col <- idx[2]
        T_result[row:(row + 1), col:(col + 1)]
      })
      T_sum <- apply(T_submatrix, 1, sum)
      R_T <- matrix(T_sum, nrow = 2) / (N * T_h)

      # Compute R_NT
      NT_de <- as.matrix(data[, .(d_e, d_lag_e)])
      R_NT <- t(NT_de) %*% NT_de / (N * T_h)

      V_TW <- solve(Q) %*% (R_N + R_T - R_NT) %*% solve(Q)
    } else {
      V_TW <- solve(Q) %*% (R_N) %*% solve(Q)
    }

    ifelse(V_TW[1, 1] < 0, 0, sqrt(V_TW[1, 1] / (N * T_h)))

  } else {
    # Compute R_N for scalar case
    N_long <- data[, c(rbind(d_e)), by = id]
    N_matrix <- matrix(N_long$V1, ncol = N)
    N_result <- N_matrix %*% t(N_matrix)
    N_sum <- sum(N_result)
    R_N <- N_sum / (N * T_h)

    # Compute variance
    if (cluster) {
      if (cluster_method == "driscoll_kraay") {
        # Driscoll and Kraay (1998) standard errors
        g <- data[, sum(d_e), by = time]$V1
        m_T <- floor(T_h^(1/4))

        # Helper function for autocovariance calculation
        calculate_Delta_j <- function(g, T_h, j) {
          valid_indices <- (j + 1):T_h
          g_t <- g[valid_indices]
          g_t_minus_j <- g[valid_indices - j]
          Delta_j <- sum(g_t * g_t_minus_j) / T_h
          return(Delta_j)
        }

        # HAC variance estimation
        Delta_0 <- calculate_Delta_j(g, T_h, j = 0)
        R_NT_DK <- Delta_0

        if (m_T >= 1) {
          w <- 1 - 1:m_T / (m_T + 1)
          for (j in 1:m_T) {
            Delta_j <- calculate_Delta_j(g, T_h, j)
            R_NT_DK <- R_NT_DK + 2 * w[j] * Delta_j
          }
        }

        R_NT_DK <- R_NT_DK / (N * T_h)
        V_TW <- solve(Q) %*% (R_N + R_NT_DK) %*% solve(Q)

      } else {
        # Two-way clustered standard errors
        # Compute R_T
        T_long <- data[, c(rbind(d_e)), by = time]
        T_matrix <- matrix(T_long$V1, ncol = T_h)
        T_result <- T_matrix %*% t(T_matrix)
        T_sum <- sum(T_result)
        R_T <- T_sum / (N * T_h)

        # Compute R_NT
        NT_de <- as.matrix(data[, .(d_e)])
        R_NT <- t(NT_de) %*% NT_de / (N * T_h)

        V_TW <- solve(Q) %*% (R_N + R_T - R_NT) %*% solve(Q)
      }

    } else {
      V_TW <- solve(Q) %*% (R_N) %*% solve(Q)
    }

    return(as.numeric(ifelse(V_TW < 0, 0, sqrt(V_TW / (N * T_h)))))
  }
}

# =============================================================================
# DGMM Standard Error Calculation
# =============================================================================

# Calculate standard errors for DGMM estimates
calculate_se_dgmm <- function(data, N, cluster = FALSE) {
  T_h <- uniqueN(data$time)

  # Get instrument and regressor matrices
  Z <- as.matrix(data[, .(iv_indep, iv_dep)])
  X <- as.matrix(data[, .(indep_diff, dep_lag_diff)])

  # Compute Q matrix
  XZ <- t(X) %*% Z
  ZZ <- solve((t(Z) %*% Z))
  ZX <- t(Z) %*% X
  Q <- XZ %*% ZZ %*% ZX / (N * T_h)

  # Compute R_N
  D_E_array <- data[, c(rbind(d_e, d_lag_e)), by = id]
  D_E <- matrix(D_E_array$V1, ncol = N)
  DEDE <- D_E %*% t(D_E)

  # Extract 2x2 submatrices
  indices <- seq(1, 2 * T_h, by = 2)
  index_grid <- expand.grid(row_start = indices, col_start = indices)
  submatrices <- apply(index_grid, 1, function(idx) {
    row <- idx[1]
    col <- idx[2]
    DEDE[row:(row + 1), col:(col + 1)]
  })
  matrices_sums <- apply(submatrices, 1, sum)

  R_N <- XZ %*% ZZ %*% matrix(matrices_sums, nrow = 2) %*% ZZ %*% ZX / (N * T_h)

  # Compute variance
  V_TW <- solve(Q) %*% (R_N) %*% solve(Q)

  ifelse(V_TW[1, 1] < 0, 0, sqrt(V_TW[1, 1] / (N * T_h)))
}
