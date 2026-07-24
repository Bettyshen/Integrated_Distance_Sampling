# Significance test among detection functions (per-panel analysis)
# Date: 2026-07-06
# Author: Betty Shen
# Goal: (1) Test whether significant pairs exist among fitted detection functions
#       (2) Identify which pairs have significant differences
#       (3) Export results to CSV aligned with figure (3 rows × 6 columns)
#
# Data structure: For each (tau, true_detection_function), compare rel_bias among
# three fitted models: half-normal, hazard-rate, negative-exponential
rm(list = ls())
# Load required libraries
library(dplyr)
library(car)   # For Levene's test (homogeneity of variance)
library(FSA)   # For Dunn's test (post-hoc for Kruskal-Wallis)

# Read results from CSV
data <- read.csv("/Volumes/T7/Chapter2/DataIntegration/Results/DataSimulation/Revision/abundance_simulation_results_compare_det.csv")
colnames(data)
# Reshape to long: one row per (pooled_replicate, tau, true_detection_function, fitted_model) with relative_bias
data_long <- bind_rows(
  data %>%
    transmute(pooled_replicate, tau, true_detection_function,
              fitted_model = "half-normal",
              relative_bias = rel_bias_half),
  data %>%
    transmute(pooled_replicate, tau, true_detection_function,
              fitted_model = "hazard-rate",
              relative_bias = rel_bias_haz),
  data %>%
    transmute(pooled_replicate, tau, true_detection_function,
              fitted_model = "negative-exponential",
              relative_bias = rel_bias_exp)
) %>%
  filter(!is.na(relative_bias)) %>%
  filter(abs(relative_bias) <= 2)  # remove outliers for all analyses

# Abbreviation helper for CSV (matches figure labels)
abbrev <- function(x) {
  case_when(
    x == "half-normal" ~ "HN",
    x == "hazard-rate" ~ "HR",
    x == "negative-exponential" ~ "NE",
    TRUE ~ x
  )
}

# Ensure factor order matches figure (rows: half-normal, hazard-rate, negative-exponential)
data_long$true_detection_function <- factor(
  data_long$true_detection_function,
  levels = c("half-normal", "hazard-rate", "negative-exponential")
)
data_long$fitted_model <- factor(
  data_long$fitted_model,
  levels = c("half-normal", "hazard-rate", "negative-exponential")
)

cat("=== Data Summary ===\n")
cat("Outliers removed: |relative_bias| > 2 excluded from all analyses.\n")
cat("Total observations (after filter):", nrow(data_long), "\n")
cat("Observations per (tau, true_detection_function):\n")
print(table(data_long$tau, data_long$true_detection_function))
cat("\n")

# ============================================================================
# Median of relative bias by scenario (aligned with figure: 18 panels)
# ============================================================================
# For each (tau, true_detection_function), median of relative_bias for each fitted_model.
median_by_scenario <- data_long %>%
  group_by(tau, true_detection_function, fitted_model) %>%
  summarise(median_rel_bias = median(relative_bias, na.rm = TRUE), .groups = "drop")

# Wide format: one row per scenario, columns median_HN, median_HR, median_NE
median_wide <- bind_rows(
  median_by_scenario %>%
    filter(fitted_model == "half-normal") %>%
    transmute(tau, true_detection_function, median_HN = median_rel_bias),
  median_by_scenario %>%
    filter(fitted_model == "hazard-rate") %>%
    transmute(tau, true_detection_function, median_HR = median_rel_bias),
  median_by_scenario %>%
    filter(fitted_model == "negative-exponential") %>%
    transmute(tau, true_detection_function, median_NE = median_rel_bias)
) %>%
  group_by(tau, true_detection_function) %>%
  summarise(
    median_HN = na.omit(median_HN)[1],
    median_HR = na.omit(median_HR)[1],
    median_NE = na.omit(median_NE)[1],
    .groups = "drop"
  ) %>%
  mutate(
    true_detection_function = factor(true_detection_function,
      levels = c("half-normal", "hazard-rate", "negative-exponential"))
  ) %>%
  arrange(true_detection_function, tau)

cat("=== MEDIAN RELATIVE BIAS BY SCENARIO ===\n")
cat("Figure layout: 3 rows (known detection function) × 6 columns (tau). HN/HR/NE = fitted model.\n\n")
print(median_wide)
cat("\n")

output_dir <- "/Volumes/T7/Chapter2/DataIntegration/Results/DataSimulation/Revision"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(median_wide, file.path(output_dir, "median_rel_bias_revision.csv"), row.names = FALSE)
cat("Median by scenario saved to: median_rel_bias_revision.csv\n\n")

# ============================================================================
# Global assumption check (full dataset — for reporting in paper)
# ============================================================================
# Single Shapiro-Wilk (residual normality) and Levene's (homogeneity of variance)
# on the full data, pooling across tau and true_detection_function.
cat("=== GLOBAL ASSUMPTION CHECK (full dataset) ===\n\n")

aov_global <- aov(relative_bias ~ fitted_model, data = data_long)
res_global <- residuals(aov_global)
n_res_global <- length(res_global)

if (n_res_global > 5000) {
  set.seed(123)
  res_global_sw <- sample(res_global, size = 5000, replace = FALSE)
  shapiro_global <- shapiro.test(res_global_sw)
  cat("Shapiro-Wilk (residual normality): n =", n_res_global, ", tested on 5000 random samples.\n")
} else {
  shapiro_global <- shapiro.test(res_global)
  cat("Shapiro-Wilk (residual normality): n =", n_res_global, ".\n")
}
cat("  W =", round(shapiro_global$statistic, 4), ", p =", format.pval(shapiro_global$p.value, digits = 3), "\n")
shapiro_global_pass <- ifelse(shapiro_global$p.value > 0.05, "Y", "N")
cat("  Pass (p > 0.05):", shapiro_global_pass, "\n\n")

levene_global <- leveneTest(relative_bias ~ fitted_model, data = data_long)
cat("Levene's test (homogeneity of variance across fitted_model groups):\n")
print(levene_global)
cat("  F =", round(levene_global$`F value`[1], 4), ", p =", format.pval(levene_global$`Pr(>F)`[1], digits = 3), "\n")
levene_global_pass <- ifelse(levene_global$`Pr(>F)`[1] > 0.05, "Y", "N")
cat("  Pass (p > 0.05):", levene_global_pass, "\n\n")

global_assump <- data.frame(
  shapiro_W = as.numeric(shapiro_global$statistic),
  shapiro_pvalue = shapiro_global$p.value,
  shapiro_pass = shapiro_global_pass,
  levene_F = levene_global$`F value`[1],
  levene_pvalue = levene_global$`Pr(>F)`[1],
  levene_pass = levene_global_pass,
  n_residuals = n_res_global,
  stringsAsFactors = FALSE
)
output_dir <- "/Volumes/T7/Chapter2/DataIntegration/Results/DataSimulation/Revision"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(global_assump, file.path(output_dir, "global_assumption_check_revision.csv"), row.names = FALSE)
cat("Global assumption check saved to:", file.path(output_dir, "global_assumption_check_revision.csv"), "\n\n")

# ============================================================================
# Per-panel analysis: Assumption checks → Kruskal-Wallis + Dunn's post-hoc
# ============================================================================
# For each panel we (1) check residuals normality (Shapiro-Wilk) and homogeneity
# of variance (Levene's test), then (2) run Kruskal-Wallis and Dunn's test.
# Figure layout: 3 rows (true detection function) × 6 columns (tau)
# Row 1: half-normal | Row 2: hazard-rate | Row 3: negative-exponential
# Col 1: tau=1 | Col 2: tau=1.5 | Col 3: tau=2 | Col 4: tau=2.5 | Col 5: tau=3 | Col 6: tau=3.5

tau_vals <- sort(unique(data_long$tau))
true_det_vals <- c("half-normal", "hazard-rate", "negative-exponential")

results_list <- list()

for (true_det in true_det_vals) {
  for (tau_val in tau_vals) {
    # Subset for this panel
    sub <- data_long %>%
      filter(tau == tau_val, true_detection_function == true_det)

    # Skip if insufficient data
    if (nrow(sub) < 6) {
      results_list[[length(results_list) + 1]] <- data.frame(
        tau = tau_val,
        true_detection_function = true_det,
        panel_row = match(true_det, true_det_vals),
        panel_col = match(tau_val, tau_vals),
        shapiro_W = NA_real_, shapiro_pvalue = NA_real_, shapiro_pass = NA_character_,
        levene_F = NA_real_, levene_pvalue = NA_real_, levene_pass = NA_character_,
        kw_chi_squared = NA,
        kw_pvalue = NA,
        kw_significant = NA,
        HN_vs_HR_pvalue = NA,
        HN_vs_HR_Z = NA_real_,
        HN_vs_HR_sig = NA,
        HN_vs_NE_pvalue = NA,
        HN_vs_NE_Z = NA_real_,
        HN_vs_NE_sig = NA,
        HR_vs_NE_pvalue = NA,
        HR_vs_NE_Z = NA_real_,
        HR_vs_NE_sig = NA,
        significant_pairs = "Insufficient data",
        n_obs = nrow(sub),
        stringsAsFactors = FALSE
      )
      next
    }

    # --- Assumption checks (before Kruskal-Wallis) ---
    # Fit ANOVA model to obtain residuals
    aov_fit <- aov(relative_bias ~ fitted_model, data = sub)
    res <- residuals(aov_fit)
    n_res <- length(res)

    # Shapiro-Wilk test for residual normality (limit 5000 as in bSim_StatTest.R)
    if (n_res > 5000) {
      set.seed(123)
      res_sw <- sample(res, size = 5000, replace = FALSE)
      shapiro_res <- shapiro.test(res_sw)
    } else {
      shapiro_res <- shapiro.test(res)
    }
    shapiro_p <- shapiro_res$p.value
    shapiro_pass <- ifelse(shapiro_p > 0.05, "Y", "N")

    # Levene's test for homogeneity of variance
    levene_res <- leveneTest(relative_bias ~ fitted_model, data = sub)
    levene_p <- levene_res$`Pr(>F)`[1]
    levene_pass <- ifelse(levene_p > 0.05, "Y", "N")

    # --- Kruskal-Wallis test ---
    kw <- kruskal.test(relative_bias ~ fitted_model, data = sub)
    kw_sig <- ifelse(kw$p.value < 0.05, "Y", "N")

    # Dunn's post-hoc
    dunn <- dunnTest(relative_bias ~ fitted_model, data = sub, method = "bonferroni")
    dunn_res <- dunn$res

    # Extract p-value and Z for each pair (Dunn format: "group1 - group2")
    get_pair <- function(g1, g2) {
      idx <- which(
        grepl(g1, dunn_res$Comparison, fixed = TRUE) & grepl(g2, dunn_res$Comparison, fixed = TRUE)
      )
      if (length(idx) > 0) {
        list(p = dunn_res$P.adj[idx[1]], z = round(dunn_res$Z[idx[1]], 2))
      } else {
        list(p = NA_real_, z = NA_real_)
      }
    }

    pair_HN_HR <- get_pair("half-normal", "hazard-rate")
    pair_HN_NE <- get_pair("half-normal", "negative-exponential")
    pair_HR_NE <- get_pair("hazard-rate", "negative-exponential")

    p_HN_HR <- pair_HN_HR$p
    p_HN_NE <- pair_HN_NE$p
    p_HR_NE <- pair_HR_NE$p

    # Significant pairs list for labeling
    sig_pairs <- c()
    if (!is.na(p_HN_HR) && p_HN_HR < 0.05) sig_pairs <- c(sig_pairs, "HN vs HR")
    if (!is.na(p_HN_NE) && p_HN_NE < 0.05) sig_pairs <- c(sig_pairs, "HN vs NE")
    if (!is.na(p_HR_NE) && p_HR_NE < 0.05) sig_pairs <- c(sig_pairs, "HR vs NE")
    sig_pairs_str <- if (length(sig_pairs) > 0) paste(sig_pairs, collapse = "; ") else "None"

    results_list[[length(results_list) + 1]] <- data.frame(
      tau = tau_val,
      true_detection_function = true_det,
      panel_row = match(true_det, true_det_vals),
      panel_col = match(tau_val, tau_vals),
      shapiro_W = shapiro_res$statistic,
      shapiro_pvalue = shapiro_p,
      shapiro_pass = shapiro_pass,
      levene_F = levene_res$`F value`[1],
      levene_pvalue = levene_p,
      levene_pass = levene_pass,
      kw_chi_squared = kw$statistic,
      kw_pvalue = kw$p.value,
      kw_significant = kw_sig,
      HN_vs_HR_pvalue = p_HN_HR,
      HN_vs_HR_Z = pair_HN_HR$z,
      HN_vs_HR_sig = ifelse(!is.na(p_HN_HR) && p_HN_HR < 0.05, "Y", "N"),
      HN_vs_NE_pvalue = p_HN_NE,
      HN_vs_NE_Z = pair_HN_NE$z,
      HN_vs_NE_sig = ifelse(!is.na(p_HN_NE) && p_HN_NE < 0.05, "Y", "N"),
      HR_vs_NE_pvalue = p_HR_NE,
      HR_vs_NE_Z = pair_HR_NE$z,
      HR_vs_NE_sig = ifelse(!is.na(p_HR_NE) && p_HR_NE < 0.05, "Y", "N"),
      significant_pairs = sig_pairs_str,
      n_obs = nrow(sub),
      stringsAsFactors = FALSE
    )
  }
}

# Combine into one data frame (rows ordered to match figure: row1 all tau, row2 all tau, row3 all tau)
results_df <- bind_rows(results_list)

# Reorder columns for clarity; drop panel_row, panel_col if not needed for labeling
results_export <- results_df %>%
  dplyr::select(tau, true_detection_function,
         shapiro_W, shapiro_pvalue, shapiro_pass,
         levene_F, levene_pvalue, levene_pass,
         kw_chi_squared, kw_pvalue, kw_significant,
         HN_vs_HR_pvalue, HN_vs_HR_Z, HN_vs_HR_sig,
         HN_vs_NE_pvalue, HN_vs_NE_Z, HN_vs_NE_sig,
         HR_vs_NE_pvalue, HR_vs_NE_Z, HR_vs_NE_sig,
         significant_pairs, n_obs)

# ============================================================================
# Export to CSV (aligned with figure: 18 rows, row order = panel order)
# ============================================================================
output_dir <- "/Volumes/T7/Chapter2/DataIntegration/Results/DataSimulation/Revision"
output_csv <- file.path(output_dir, "stat_test_results_significance_test.csv")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(results_export, file = output_csv, row.names = FALSE)

# ============================================================================
# Export pairwise Dunn results (long format for tables/figures)
# Columns: tau value, true detection function, comparison pairs, Z score, p value
# ============================================================================
format_pvalue_export <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("p<0.001")
  sprintf("%.2f", p)
}

true_det_label <- function(x) {
  case_when(
    x == "half-normal" ~ "Half-normal",
    x == "hazard-rate" ~ "Hazard-rate",
    x == "negative-exponential" ~ "Negative-exponential",
    TRUE ~ as.character(x)
  )
}

pairwise_export <- bind_rows(
  results_export %>%
    transmute(
      `tau value` = tau,
      `true detection function` = true_det_label(true_detection_function),
      `detection function comparison pairs` = "HN vs HR",
      `Z score` = round(HN_vs_HR_Z, 2),
      `p value` = vapply(HN_vs_HR_pvalue, format_pvalue_export, character(1))
    ),
  results_export %>%
    transmute(
      `tau value` = tau,
      `true detection function` = true_det_label(true_detection_function),
      `detection function comparison pairs` = "HN vs NE",
      `Z score` = round(HN_vs_NE_Z, 2),
      `p value` = vapply(HN_vs_NE_pvalue, format_pvalue_export, character(1))
    ),
  results_export %>%
    transmute(
      `tau value` = tau,
      `true detection function` = true_det_label(true_detection_function),
      `detection function comparison pairs` = "HR vs NE",
      `Z score` = round(HR_vs_NE_Z, 2),
      `p value` = vapply(HR_vs_NE_pvalue, format_pvalue_export, character(1))
    )
) %>%
  arrange(
    factor(`true detection function`,
           levels = c("Half-normal", "Hazard-rate", "Negative-exponential")),
    `tau value`,
    factor(`detection function comparison pairs`,
           levels = c("HN vs HR", "HN vs NE", "HR vs NE"))
  )

pairwise_csv <- file.path(output_dir, "stat_test_pairwise_dunn_revision.csv")
write.csv(pairwise_export, file = pairwise_csv, row.names = FALSE)

cat("=== STATISTICAL TEST RESULTS (per panel) ===\n\n")
cat("Figure layout: 3 rows (known detection function) × 6 columns (tau)\n")
cat("Each row below = one panel, ordered top-to-bottom, left-to-right.\n\n")
print(results_export)
cat("\n")
cat("Exported to:", output_csv, "\n")
cat("Pairwise Dunn export saved to:", pairwise_csv, "\n\n")


# ============================================================================
# Perceptability summary (by tau × true_detection_function)
#   Note: use original wide data since it contains `perceptability`
# ============================================================================
perceptability_summary <- data %>%
  group_by(tau, true_detection_function) %>%
  summarise(
    mean_perceptability   = mean(perceptability, na.rm = TRUE),
    median_perceptability = median(perceptability, na.rm = TRUE),
    q25_perceptability    = quantile(perceptability, probs = 0.25, na.rm = TRUE),
    q75_perceptability    = quantile(perceptability, probs = 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(perceptability_summary, file.path(output_dir, "perceptability_summary_revision.csv"), row.names = FALSE)
cat("Perceptability summary saved to:", file.path(output_dir, "perceptability_summary_revision.csv"), "\n\n")

# ============================================================================
# Console summary (assumption checks + test results)
# ============================================================================
cat("=== ASSUMPTION CHECK SUMMARY (per panel) ===\n")
cat("Shapiro-Wilk (normality of residuals): pass = p > 0.05. Levene (homogeneity of variance): pass = p > 0.05.\n\n")
for (i in seq_len(nrow(results_export))) {
  r <- results_export[i, ]
  cat(sprintf("Panel [%s, tau=%.3g]: Shapiro p=%.2e (%s) | Levene p=%.4f (%s)\n",
              r$true_detection_function, r$tau,
              r$shapiro_pvalue, r$shapiro_pass,
              r$levene_pvalue, r$levene_pass))
}
cat("\n=== KRUSKAL-WALLIS & PAIRWISE SUMMARY ===\n\n")
for (i in seq_len(nrow(results_export))) {
  r <- results_export[i, ]
  cat(sprintf("Panel [%s, tau=%.3g]: KW p=%.4f (%s) | Sig pairs: %s\n",
              r$true_detection_function, r$tau,
              r$kw_pvalue, r$kw_significant,
              r$significant_pairs))
}
cat("\n=== ANALYSIS COMPLETE ===\n")
