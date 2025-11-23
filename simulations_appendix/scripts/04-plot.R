# =============================================================================
# 04-plot.r
# Purpose: Create and save visualizations for simulation results
# =============================================================================

# =============================================================================
# Plotting Theme
# =============================================================================

# Define a reusable plotting theme
# Creates a consistent theme for ggplot objects with customized aesthetics
plot_with_theme <- function(base_family = "Times New Roman", base_size = 13) {
  theme_bw(base_family = base_family, base_size = base_size) +
    theme(
      # Facet strip settings
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 18),
      strip.text.y = element_text(size = 18, angle = 90),
      # Axis settings
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      # Panel settings
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "white", linewidth = 0.7),
      panel.spacing = unit(0.35, "in"),
      # Legend settings
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 14, face = "plain"),
      legend.key.spacing = unit(0.4, "in"),
      legend.key.width = unit(0.4, "in"),
      # Other settings
      aspect.ratio = 1
    )
}

# =============================================================================
# Facet Functions
# =============================================================================

# Define facet grid layout with AR coefficient rows and N, T columns
facet_with_labeller <- function() {
  facet_grid(
    AR_Coefficient ~ N + T,
    labeller = label_bquote(
      cols = N == .(N) ~ "," ~ T == .(T),
      rows = rho == .(AR_Coefficient)
    ),
    switch = "y",
    axes = "all"
  )
}

# Define facet grid layout for contrast plots (N, T columns only)
facet_with_labeller_contrast <- function() {
  facet_grid(
    cols = vars(N, T),
    labeller = label_bquote(
      cols = N == .(N) ~ "," ~ T == .(T)
    ),
    axes = "all"
  )
}

# =============================================================================
# Plot Saving Functions
# =============================================================================

# Save plot with embedded fonts for cross-system compatibility
save_plot_with_theme <- function(plot, filename, width = 12, height = 15, embed_fonts = TRUE) {
  ggsave(filename, plot, width = width, height = height, units = "in")
  if (embed_fonts) {
    embed_fonts(filename)
  }
}

# =============================================================================
# Plot Creation Functions
# =============================================================================

# Create customized plot with consistent styling
# Handles data processing, aesthetics, and optional reference lines
create_plot <- function(data, x, y, estimator, facet_func, scale_y = NULL, filename, y_intercept = NULL, hline_type = "dashed") {
  plot <- data |>
    ggplot(
      aes(x = {{ x }}, y = {{ y }}, color = {{ estimator }}, shape = {{ estimator }})
    ) +
    geom_line(linewidth = 0.4, alpha = 0.4) +
    geom_point(size = 1, alpha = 0.4) +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    # Estimator shape mapping
    scale_shape_manual(values = c(
      "FE" = 1, "SPJ" = 2, "true IRF" = 0, "true value" = 0,
      "DGMM" = 5, "IVXJ" = 5, "DB" = 5, "LE" = 4
    )) +
    # Estimator color mapping
    scale_color_manual(values = c(
      "IVXJ" = "black", "DGMM" = "black", "FE" = "blue",
      "SPJ" = "firebrick3", "true IRF" = "#009900",
      "true value" = "#009900", "DB" = "black", "LE" = "orange"
    )) +
    facet_func() +
    plot_with_theme() +
    guides(
      color = guide_legend(override.aes = list(size = 3, alpha = 1)),
      linetype = guide_legend(override.aes = list(size = 2, alpha = 1))
    )

  # Add optional horizontal reference line
  if (!is.null(y_intercept)) {
    plot <- plot + geom_hline(yintercept = y_intercept, linetype = hline_type)
  }

  # Apply y-axis scale if specified
  if (!is.null(scale_y)){
    plot <- plot +
      scale_y_continuous(limits = scale_y$limits, breaks = scale_y$breaks)
  }

  save_plot_with_theme(plot, filename)
}

# =============================================================================
# Main Plotting Function
# =============================================================================

# Generate and save multiple plots for simulation results
plot_simulation_results <- function(data, limits_vs, breaks_vs, case, vars_num = 2) {

  if (vars_num == 2) {
    # Create IRF plot
    create_plot(
      data, Horizon, beta, estimator,
      facet_func = facet_with_labeller,
      scale_y = list(limits = limits_vs, breaks = breaks_vs),
      filename = paste0("results/plots/", case, "_IRF_Estimates_vs_True_IRFs.pdf"),
      y_intercept = 0
    )

    # Create RMSE plot
    data <- data |>
      filter(estimator != "true IRF")
    create_plot(
      data, Horizon, rmse, estimator,
      facet_func = facet_with_labeller,
      scale_y = list(limits = c(0, 0.1002), breaks = c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10)),
      filename = paste0("results/plots/", case, "_RMSE.pdf")
    )

    # Create coverage probability plot
    create_plot(
      data, Horizon, coverage, estimator,
      facet_func = facet_with_labeller,
      scale_y = list(limits = c(0, 1), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
      filename = paste0("results/plots/", case, "_Coverage_Probabilities.pdf"),
      y_intercept = 0.95
    )

    # Create confidence interval lengths plot
    create_plot(
      data, Horizon, confidence_interval_length, estimator,
      facet_func = facet_with_labeller,
      scale_y = list(limits = c(0, 0.23), breaks = c(0.00, 0.05, 0.10, 0.15, 0.20)),
      filename = paste0("results/plots/", case, "_confidence_interval_length.pdf")
    )

  } else {
    # Create IRF plot (contrast layout)
    create_plot(
      data, Horizon, beta, estimator,
      facet_func = facet_with_labeller_contrast,
      filename = paste0("results/plots/", case, "_IRF_Estimates_vs_True_IRFs.pdf"),
      y_intercept = 0
    )

    # Create RMSE plot (contrast layout)
    data2 <- data |>
      filter(estimator != "true IRF")
    create_plot(
      data2, Horizon, rmse, estimator,
      facet_func = facet_with_labeller_contrast,
      scale_y = list(limits = c(0, 0.1002), breaks = c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10)),
      filename = paste0("results/plots/", case, "_RMSE.pdf")
    )

    # Create contrast plot
    data_contrast <- data |>
      mutate(estimator = ifelse(estimator == "true IRF", "true value", estimator))
    create_plot(
      data_contrast, Horizon, contrast_YZ, estimator,
      facet_func = facet_with_labeller_contrast,
      filename = paste0("results/plots/", case, "_Contrast.pdf")
    )
  }
}
