# Abundance Comparison Simulation (equal availability)
# Goal: Compare the correct detection function vs. incorrect misspecified detection functions results in density estimates.
# Date: 2026-07-04
# Author: Betty Shen
# Correct vs incorrect detection functions in distsamp.
rm(list = ls())
# Load required libraries
library(bSims)
library(detect)
library(unmarked)
library(dplyr)

# Set simulation parameters
n_sites_per_pool <- 30     # point counts pooled into one distsamp fit (rows in y matrix)
n_pooled_replicates <- 500  # independent pooled fits per tau × true detection function (iterations)
tau_values <- c(1, 1.5, 2, 2.5, 3, 3.5)
phi <- 1  # recorded in output only; ignored when initial_location = TRUE
b <- 3 # hazard-rate shape parameter (shoulder width)
availability <- 1

# Landscape and detection parameters
extent <- 10  # 1 km x 1 km landscape
density <- 1  # 1 individual/hectare = 100 individuals total
duration <- 5  # survey duration label (events fixed at t = 0 with initial_location)
tbr <- c(5)
rbr <- seq(0.5, 4, by = 0.5)           # distance bins (bSims rint units)
dist_breaks <- c(0, rbr)             # Ken-style breaks with unitsIn = "km"

half_normal_fun <- function(d, tau) exp(-d^2 / (2 * tau^2))
exponential_fun <- function(d, tau) exp(-d / tau)
hazard_fun <- function(d, tau) 1 - exp(-(d / tau)^-b)

true_det_specs <- list(
  list(label = "half-normal", dist_fun = half_normal_fun),
  list(label = "hazard-rate", dist_fun = hazard_fun),
  list(label = "negative-exponential", dist_fun = exponential_fun)
)

# Initialize results storage
results <- data.frame(
  pooled_replicate = integer(),
  true_detection_function = character(),
  n_sites = integer(),
  true_abundance = numeric(),
  est_half = numeric(),
  est_haz = numeric(),
  est_exp = numeric(),
  diff_true_half = numeric(),
  diff_true_haz = numeric(),
  diff_true_exp = numeric(),
  diff_half_haz = numeric(),
  diff_half_exp = numeric(),
  diff_haz_exp = numeric(),
  rel_bias_half = numeric(),
  rel_bias_haz = numeric(),
  rel_bias_exp = numeric(),
  availability = numeric(),
  perceptability = numeric(),
  tau = numeric(),
  b = numeric(),
  phi = numeric(),
  stringsAsFactors = FALSE
)

# One independent point count in bSims → distance-bin counts + site-level summaries
run_point_count <- function(dist_fun, tau) {
  l <- bsims_init(extent = extent, road = 0, edge = 0)
  p <- bsims_populate(l, density = density) # density = 1 individual/hectare = 100 individuals total
  e <- bsims_animate(p, initial_location = TRUE, duration = duration) # duration = 5 minutes; initial_location = TRUE makes all birds equally available (available = 1)
  true_abundance <- get_abundance(p) # true number of individuals in landscape
  d <- bsims_detect(e, dist_fun = dist_fun, tau = tau, event_type = "vocal") # only consider singing individuals
  x <- bsims_transcribe(d, tint = tbr, rint = rbr)
  counts <- rowSums(get_table(x, type = "removal"))
  perceptability <- function(dist_fun, tau, w = 4) {
  (2 / w^2) *
    integrate(
      function(x) x * dist_fun(x, tau),
      lower = 0,
      upper = w
    )$value
}
  list(counts = counts, true_abundance = true_abundance, perceptability = perceptability)
}

# Pooled fit in unmarked distance sampling: y_matrix has one row per point count (site)
fit_pooled_keyfuns <- function(y_matrix, mean_abundance, mean_perceptability,
                               pooled_replicate, tau, true_detection_label) {
  umf <- unmarkedFrameDS(
    y = y_matrix,
    dist.breaks = dist_breaks,
    survey = "point",
    unitsIn = "km"
  )

  fit_half <- distsamp(~1 ~1, data = umf, keyfun = "halfnorm", unitsOut = "km")
  fit_haz <- distsamp(~1 ~1, data = umf, keyfun = "hazard", unitsOut = "km")
  fit_exp <- distsamp(~1 ~1, data = umf, keyfun = "exp", unitsOut = "km")

# backTransform is from unmarked package and convert individuals units to per km^2
  est_half <- backTransform(fit_half, type = "state")@estimate * 100
  est_haz <- backTransform(fit_haz, type = "state")@estimate * 100
  est_exp <- backTransform(fit_exp, type = "state")@estimate * 100

  data.frame(
    pooled_replicate = pooled_replicate,
    true_detection_function = true_detection_label,
    n_sites = nrow(y_matrix),
    true_abundance = mean_abundance,
    est_half = est_half,
    est_haz = est_haz,
    est_exp = est_exp,
    diff_true_half = mean_abundance - est_half,
    diff_true_haz = mean_abundance - est_haz,
    diff_true_exp = mean_abundance - est_exp,
    diff_half_haz = est_half - est_haz,
    diff_half_exp = est_half - est_exp,
    diff_haz_exp = est_haz - est_exp,
    rel_bias_half = (est_half - mean_abundance) / mean_abundance,
    rel_bias_haz = (est_haz - mean_abundance) / mean_abundance,
    rel_bias_exp = (est_exp - mean_abundance) / mean_abundance,
    availability = availability,
    perceptability = mean_perceptability,
    tau = tau,
    b = b,
    phi = phi,
    stringsAsFactors = FALSE
  )
}

na_result_row <- function(pooled_replicate, tau, mean_abundance, true_detection_label) {
  data.frame(
    pooled_replicate = pooled_replicate,
    true_detection_function = true_detection_label,
    n_sites = n_sites_per_pool,
    true_abundance = mean_abundance,
    est_half = NA, est_haz = NA, est_exp = NA,
    diff_true_half = NA, diff_true_haz = NA, diff_true_exp = NA,
    diff_half_haz = NA, diff_half_exp = NA, diff_haz_exp = NA,
    rel_bias_half = NA, rel_bias_haz = NA, rel_bias_exp = NA,
    availability = availability,
    perceptability = NA,
    tau = tau, b = b, phi = phi,
    stringsAsFactors = FALSE
  )
}

# Expected number of rows in results dataframe
expected_rows <- length(tau_values) * length(true_det_specs) * n_pooled_replicates
cat("Sites per pooled fit:", n_sites_per_pool, "\n")
cat("Iterations per scenario:", n_pooled_replicates, "\n")
cat("Tau values:", paste(tau_values, collapse = ", "), "\n")
cat("Expected result rows:", expected_rows, "\n\n")

# Run simulations for each tau value and true detection function
#set.seed(123)

for (tau in tau_values) {
  cat("\n=== Running simulations for tau =", tau, "===\n")

  for (spec in true_det_specs) {
    cat("  true detection:", spec$label, "\n")

    for (rep in seq_len(n_pooled_replicates)) {
      site_runs <- replicate(n_sites_per_pool, run_point_count(spec$dist_fun, tau),
                             simplify = FALSE)
      y_matrix <- do.call(rbind, lapply(site_runs, `[[`, "counts"))
      mean_abundance <- mean(vapply(site_runs, `[[`, numeric(1), "true_abundance"))
      mean_perceptability <- mean(vapply(site_runs, `[[`, numeric(1), "perceptability"))

      tryCatch({
        results <- rbind(
          results,
          fit_pooled_keyfuns(y_matrix, mean_abundance, mean_perceptability,
                             rep, tau, spec$label)
        )
      }, error = function(err) {
        results <<- rbind(
          results,
          na_result_row(rep, tau, mean_abundance, spec$label)
        )
      })
    }
  }

  cat("Completed tau =", tau, "\n")
}

cat("\n=== All simulations complete! ===\n")

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total rows:", nrow(results), "\n")
cat("Expected rows:", expected_rows, "\n")

cat("\nRows by tau:\n")
print(table(results$tau))

cat("\nRows by true detection function:\n")
print(table(results$true_detection_function))

cat("\nMean mismatch metrics by tau and true detection function:\n")
summary_table <- aggregate(
  cbind(
    mean_diff_true_half = diff_true_half,
    mean_diff_true_haz = diff_true_haz,
    mean_diff_true_exp = diff_true_exp,
    mean_rel_bias_half = rel_bias_half,
    mean_rel_bias_haz = rel_bias_haz,
    mean_rel_bias_exp = rel_bias_exp
  ) ~ tau + true_detection_function,
  data = results,
  FUN = function(x) mean(x, na.rm = TRUE)
)
print(summary_table)

# Export results to CSV
output_file <- "/Volumes/T7/Results/DataSimulation/Revision/abundance_simulation_results_compare_det.csv"
write.csv(results, file = output_file, row.names = FALSE)
cat("\nResults exported to:", output_file, "\n")

cat("\nFirst 10 rows of results:\n")
print(head(results, 10))

# Perceptability summary
perceptability_summary <- results %>%
  group_by(tau, true_detection_function) %>%
  summarise(
    mean_perceptability = round(mean(perceptability, na.rm = TRUE), 2),
    median_perceptability = round(median(perceptability, na.rm = TRUE), 2),
    q25_perceptability = round(quantile(perceptability, probs = 0.25, na.rm = TRUE), 2),
    q75_perceptability = round(quantile(perceptability, probs = 0.75, na.rm = TRUE), 2),
    .groups = "drop"
  )
print(perceptability_summary)

perceptability_file <- "/Volumes/T7/DataSimulation/Revision/perceptability_summary_compare_det.csv"
write.csv(perceptability_summary, file = perceptability_file, row.names = FALSE)
cat("\nPerceptability summary exported to:", perceptability_file, "\n")

# Relative bias medians
rel_bias_medians <- results %>%
  group_by(tau, true_detection_function) %>%
  summarise(
    median_rel_bias_half = round(median(rel_bias_half, na.rm = TRUE), 2),
    median_rel_bias_haz = round(median(rel_bias_haz, na.rm = TRUE), 2),
    median_rel_bias_exp = round(median(rel_bias_exp, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    half_label = if_else(true_detection_function == "half-normal", "correct", "non-correct"),
    haz_label = if_else(true_detection_function == "hazard-rate", "correct", "non-correct"),
    exp_label = if_else(true_detection_function == "negative-exponential", "correct", "non-correct")
  )

rel_bias_long <- bind_rows(
  rel_bias_medians %>%
    transmute(tau, true_detection_function,
              fitted_model = "half-normal",
              median_rel_bias = median_rel_bias_half,
              correct_or_not = half_label),
  rel_bias_medians %>%
    transmute(tau, true_detection_function,
              fitted_model = "hazard-rate",
              median_rel_bias = median_rel_bias_haz,
              correct_or_not = haz_label),
  rel_bias_medians %>%
    transmute(tau, true_detection_function,
              fitted_model = "negative-exponential",
              median_rel_bias = median_rel_bias_exp,
              correct_or_not = exp_label)
)

cat("\n=== RELATIVE BIAS MEDIANS ===\n")
print(rel_bias_long)

rel_bias_file <- "/Volumes/T7/DataIntegration/Results/DataSimulation/Revision/rel_bias_medians_compare_det.csv"
write.csv(rel_bias_long, file = rel_bias_file, row.names = FALSE)
cat("\nRelative bias medians exported to:", rel_bias_file, "\n")
