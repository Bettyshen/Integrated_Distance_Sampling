# Abundance Comparison Simulation - Violin plots
# 12 panels: 3 rows (true detection function) × 4 columns (tau)
# Row 1: half-normal is correct; Row 2: hazard-rate is correct; Row 3: negative-exponential is correct
# Date: 2026-03-17

# Load required libraries
library(ggplot2)
library(dplyr)
library(grid)
library(gtable)

# Read output CSV files
results <- read.csv("/Volumes/T7/Chapter2/DataIntegration/Results/DataSimulation/PointCounts25/density_simulation_results_25counts.csv")
perceptibility <- read.csv("/Volumes/T7/Chapter2/DataIntegration/Results/DataSimulation/PointCounts25/perceptability_summary_25counts.csv")

# Tau values used in the simulation (must match abundance_simulation_*_para.R)
tau_values <- c(1.5, 2, 2.5, 3)

# Reshape results to long: one row per (tau, true_detection_function, simulation, fitted_model) with relative_bias
results_long <- bind_rows(
  results %>%
    transmute(tau, true_detection_function, simulation,
              fitted_model = "half-normal",
              relative_bias = rel_bias_half),
  results %>%
    transmute(tau, true_detection_function, simulation,
              fitted_model = "hazard-rate",
              relative_bias = rel_bias_haz),
  results %>%
    transmute(tau, true_detection_function, simulation,
              fitted_model = "negative-exponential",
              relative_bias = rel_bias_exp)
) %>%
  mutate(correct_or_not = (fitted_model == true_detection_function))

# Exclude relative_bias values with absolute value > 2 for clearer violin display
results_long <- results_long %>%
  filter(is.na(relative_bias) | abs(relative_bias) <= 2)

# Factor order: rows = true detection function (half-normal, hazard-rate, negative-exponential)
# Columns = tau 1, 2, 3, 4
results_long$tau <- factor(results_long$tau, levels = tau_values)
results_long$true_detection_function <- factor(
  results_long$true_detection_function,
  levels = c("half-normal", "hazard-rate", "negative-exponential")
)
results_long$fitted_model <- factor(
  results_long$fitted_model,
  levels = c("half-normal", "hazard-rate", "negative-exponential")
)

# Perceptibility annotation: one value per panel (tau × true_detection_function)
perceptibility$tau <- factor(perceptibility$tau, levels = tau_values)
perceptibility$true_detection_function <- factor(
  perceptibility$true_detection_function,
  levels = c("half-normal", "hazard-rate", "negative-exponential")
)
perceptibility_annot <- perceptibility %>%
  mutate(y_pos = 2.3)

# 12-panel plot: facet_grid(true_detection_function ~ tau), 3 rows × 4 columns
p <- results_long %>%
  ggplot(aes(x = fitted_model, y = relative_bias, fill = fitted_model)) +
  geom_violin(width = 0.85, scale = "width") +
  geom_boxplot(width = 0.12, color = "white", alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_text(
    data = perceptibility_annot,
    aes(x = 2, y = y_pos, label = sprintf("(%.2f)", median_perceptability)),
    inherit.aes = FALSE,
    size = 8,
    hjust = 0.5,
    vjust = 0
  ) +
  scale_fill_manual(
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
  scale_x_discrete(labels = c(
    "half-normal" = "HN",
    "hazard-rate" = "HR",
    "negative-exponential" = "NE"
  )) +
  scale_y_continuous(limits = c(NA, 2.5)) +
  facet_grid(true_detection_function ~ tau,
             switch = "y",
             labeller = labeller(
               tau = function(x) paste("tau =", x),
               true_detection_function = function(x) {
                 x <- gsub("-", " ", x)
                 case_when(
                   x == "half normal" ~ "(A) Half-normal",
                   x == "hazard rate" ~ "(B) Hazard rate",
                   x == "negative exponential" ~ "(C) Negative exponential",
                   TRUE ~ x
                 )
               }
             )) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 25, face = "bold"),
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text = element_text(size = 22, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20, margin = margin(b = 8, t = 8, l = 12, unit = "pt")),
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.key.spacing.y = unit(0.3, "cm"),
    axis.title = element_text(size = 25),
    axis.title.x = element_text(margin = margin(t = 20, unit = "pt"), face = "bold"),
    axis.title.y = element_text(margin = margin(r = 20, unit = "pt"), face = "bold"),
    axis.text = element_text(size = 17),
    axis.text.x = element_text(size = 20),
    panel.spacing = unit(1, "lines"),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.075, 0.99)  # nudge "Known detection function" right from topleft (x, y in npc)
  ) +
  labs(fill = "Fitted detection function",
       tag = "") +
  xlab("Fitted detection function") +
  ylab(NULL)

# Add a dedicated column between facet strip and panels, then place "Relative bias" in it
gt <- ggplotGrob(p)
strip_idx <- which(grepl("strip-l", gt$layout$name, fixed = TRUE))
panel_idx <- which(grepl("^panel", gt$layout$name))  # "panel-1-1", "panel-2-1", etc.
if (length(strip_idx) > 0 && length(panel_idx) > 0) {
  strip_col <- max(gt$layout$r[strip_idx])
  t_panel <- min(gt$layout$t[panel_idx])
  b_panel <- max(gt$layout$b[panel_idx])
  # Insert a column of space (1.5 cm) between strip and panels
  gt <- gtable_add_cols(gt, widths = unit(1.5, "cm"), pos = strip_col)
  rel_bias_grob <- textGrob("Relative bias", rot = 90,
                            gp = gpar(fontsize = 25, fontface = "bold"))
  gt <- gtable_add_grob(gt, rel_bias_grob,
                        t = t_panel, b = b_panel,
                        l = strip_col + 1L, r = strip_col + 1L,
                        name = "ylab_rel_bias", clip = "off")
  p_labeled <- gt
} else {
  p_labeled <- p
}

# Save plot (3 rows × 4 columns = 12 panels)
ggsave("/Volumes/T7/Chapter2/DataIntegration/Results/DataSimulation/PointCounts25/Plot/relative_bias_violin_plot_tau_panels_publish_25counts.png",
       plot = p_labeled,
       width = 25,
       height = 16,
       dpi = 300)
