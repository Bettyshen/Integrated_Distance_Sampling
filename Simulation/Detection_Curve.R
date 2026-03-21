library(ggplot2)

# Distance (in metres) out to truncation
d <- seq(0, 400, by = 1)

# Tau values correspond to the simulation settings:
#   tau_values <- c(0.1, 0.25, 0.5, 1)
# and are converted to metres via tau_m = tau_value * 100
tau_units <-  c(1.5, 2, 2.5, 3)
tau_m_vals <- tau_units * 100

# Hazard-rate shape parameters (controls shoulder)
b_vals <- c(3)

# Build a data frame of detection curves for:
#  - half-normal
#  - negative-exponential
#  - hazard-rate with varying b
curves <- list()

idx <- 1
for (tau_m in tau_m_vals) {
  for (dist_fun in c("half-normal", "negative-exponential")) {
    if (dist_fun == "half-normal") {
      g <- exp(- (d^2) / (2 * tau_m^2))
    } else if (dist_fun == "negative-exponential") {
      g <- exp(-(d / tau_m))
    }

    curves[[idx]] <- data.frame(
      distance = d,
      g = g,
      tau_unit = tau_m / 100,
      tau_m = tau_m,
      kernel = dist_fun,
      b = NA_real_
    )
    idx <- idx + 1
  }

  for (b in b_vals) {
    g <- 1 - exp(-(d / tau_m)^(-b))
    curves[[idx]] <- data.frame(
      distance = d,
      g = g,
      tau_unit = tau_m / 100,
      tau_m = tau_m,
      kernel = "hazard-rate",
      b = b
    )
    idx <- idx + 1
  }
}

curve_df <- do.call(rbind, curves)

curve_df$tau_f <- factor(
  curve_df$tau_unit,
  levels = tau_units,
  labels = paste0("tau = ", tau_units)
)

output_dir <- "/Volumes/T7/Chapter2/DataIntegration/Results/DataSimulation/PointCounts25"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
png_file <- file.path(output_dir, "detection_curves_tau_b.png")

png(png_file, width = 2800, height = 1600, res = 300)

print(
  ggplot(curve_df, aes(x = distance, y = g, colour = kernel, linetype = kernel)) +
    geom_line(size = 0.8) +
    facet_wrap(~ tau_f, nrow = 2) +
    scale_x_continuous("Distance (m)", limits = c(0, 400), breaks = seq(0, 400, 100)) +
    scale_y_continuous("Perceptibility", limits = c(0, 1)) +
    scale_colour_manual(
      values = c(
        "half-normal" = "#E69F00",
        "hazard-rate" = "#56B4E9",
        "negative-exponential" = "#009E73"
      ),
      labels = c(
        "half-normal" = "Half-normal (HN)",
        "hazard-rate" = "Hazard rate (HR)",
        "negative-exponential" = "Negative exponential (NE)"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "half-normal" = "solid",
        "hazard-rate" = "solid",
        "negative-exponential" = "solid"
      )
    ) +
    guides(linetype = "none") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14, face = "plain"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 14, face = "bold")
    ) 
)

dev.off()

cat("Detection curves saved to:", png_file, "\n")
