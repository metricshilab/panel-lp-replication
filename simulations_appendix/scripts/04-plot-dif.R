# =============================================================================
# 04-plot-dif.r
# Purpose: Create and save visualizations for difference-based simulation results
# =============================================================================

# Load required libraries
library(gridExtra)
library(patchwork)

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
      aspect.ratio = 1,
      plot.caption = element_text(size = 14, hjust = 0.5)
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
    ) +
    theme(
      strip.text.y = element_blank(),
      strip.text.x = element_blank(),
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    labs(title = filename)

  # Add optional horizontal reference line
  if (!is.null(y_intercept)) {
    plot <- plot + geom_hline(yintercept = y_intercept, linetype = hline_type)
  }

  # Apply y-axis scale if specified
  if (!is.null(scale_y)){
    plot <- plot +
      scale_y_continuous(breaks = scale_y$breaks) +
      coord_cartesian(ylim = scale_y$limits)
  }

  return(plot)
}

# =============================================================================
# Main Plotting Function
# =============================================================================

# Generate and save multiple plots for simulation results
plot_simulation_results <- function(data, limits_vs, breaks_vs, case, vars_num = 2) {

  if (vars_num == 2) {
    # Create IRF plot
    plot1 <- create_plot(
      data, Horizon, beta, estimator,
      facet_func = facet_with_labeller,
      scale_y = list(limits = limits_vs, breaks = breaks_vs),
      filename = "Averaged Estimated IRFs",
      y_intercept = 0
    )

    # Create RMSE plot
    data <- data |>
      filter(estimator != "true IRF")
    plot2 <- create_plot(
      data, Horizon, rmse, estimator,
      facet_func = facet_with_labeller,
      scale_y = list(limits = c(0, 0.12), breaks = c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12)),
      filename = "RMSEs"
    ) + theme(legend.position = "none")

    # Create coverage probability plot
    plot3 <- create_plot(
      data, Horizon, coverage, estimator,
      facet_func = facet_with_labeller,
      scale_y = list(limits = c(0, 1), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
      filename = "Coverage Probabilities",
      y_intercept = 0.95
    ) + theme(legend.position = "none")

    # Combine plots with shared legend at bottom
    combined_plot <- (plot1 | plot2 | plot3) +
      plot_layout(guides = "collect") +
      plot_annotation(
        theme = theme(legend.position = "bottom")
      )

    # Save combined plot
    ggsave("results/plots/dif.pdf", combined_plot, width = 12, height = 6)
  }
}
