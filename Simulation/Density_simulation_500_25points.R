# Abundance Comparison Simulation
# Date: 2026-03-16
# 25 points counts per simulation
# 500 simulations per tau (4 tau values), 3 true detection functions → 6000 datasets
# correct vs incorrect detection functions in distsamp.

library(bSims)
library(detect)
library(unmarked)
library(parallel)

# Simulation parameters
n_simulations <- 500
tau_values <- c(1.5, 2, 2.5, 3)
n_sites <- 25

phi <- 0.3
b <- 3
availability <- 0.4

# Landscape parameters
extent <- 10
density <- 1
duration <- 5

tbr <- seq(from = 1, to = 5, by = 1)
rbr <- seq(from = 0.5, to = 4, by = 0.5)

dist.breaks.meters <- c(0, rbr * 100)

# Initialize results storage
results <- data.frame(
  simulation = integer(),
  true_detection_function = character(),
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

# Prepare y matrix for unmarked
prepare_for_unmarked <- function(transcribe_list) {

  y_list <- lapply(transcribe_list, function(x) {

    tab <- get_table(x, type = "removal")

    rowSums(tab)

  })

  y_matrix <- do.call(rbind, y_list)

  return(y_matrix)

}

# Calculate perceptibility
# detections_list: list of bsims_detections objects (one per site/point count)
# true_abundance: integer (true number of individuals in landscape)
# verbose: logical; if TRUE, prints first detection column names for debugging

# Per-point-count average perceptibility
calculate_perceptability <- function(detections_list,
                                      true_abundance,
                                      verbose = FALSE) {
  if (length(detections_list) == 0 || true_abundance <= 0) return(0)

  per_site <- numeric(length(detections_list))

  for (s in seq_along(detections_list)) {
    x <- detections_list[[s]]

    det <- tryCatch(get_detections(x), error = function(e) NULL)
    if (is.null(det) || nrow(det) == 0) {
      per_site[s] <- 0
      next
    }

    if (verbose && s == 1) {
      message("detect columns (first site): ", paste(names(det), collapse = ", "))
    }
    # Match original approach: use number of detected rows as a perceptibility proxy
    per_site[s] <- nrow(det) / true_abundance
  }

  # average over point counts, clamp to [0,1]
  p <- mean(per_site)
  p
}

# Model fitting function
fit_all_keyfuns <- function(y,
                            true_abundance,
                            sim,
                            tau,
                            true_detection_label,
                            siteCovs,
                            perceptability_value) {

  umf <- unmarkedFrameDS(
    y = y,
    siteCovs = siteCovs,
    dist.breaks = dist.breaks.meters,
    survey = "point",
    unitsIn = "m"
  )

  fit_half <- distsamp(~1 ~1, data = umf, keyfun = "halfnorm")
  fit_haz  <- distsamp(~1 ~1, data = umf, keyfun = "hazard")
  fit_exp  <- distsamp(~1 ~1, data = umf, keyfun = "exp")

  est_half <- predict(fit_half, type = "state")[1,1] * 100
  est_haz  <- predict(fit_haz, type = "state")[1,1] * 100
  est_exp  <- predict(fit_exp, type = "state")[1,1] * 100

  data.frame(
    simulation = sim,
    true_detection_function = true_detection_label,
    true_abundance = true_abundance,

    est_half = est_half,
    est_haz = est_haz,
    est_exp = est_exp,

    diff_true_half = true_abundance - est_half,
    diff_true_haz  = true_abundance - est_haz,
    diff_true_exp  = true_abundance - est_exp,

    diff_half_haz = est_half - est_haz,
    diff_half_exp = est_half - est_exp,
    diff_haz_exp  = est_haz - est_exp,

    rel_bias_half = ifelse(true_abundance > 0,
                           (est_half - true_abundance)/true_abundance, NA),

    rel_bias_haz = ifelse(true_abundance > 0,
                          (est_haz - true_abundance)/true_abundance, NA),

    rel_bias_exp = ifelse(true_abundance > 0,
                          (est_exp - true_abundance)/true_abundance, NA),

    availability = availability,
    perceptability = perceptability_value,
    tau = tau,
    b = b,
    phi = phi,
    stringsAsFactors = FALSE
  )

}

# NA row helper
na_result_row <- function(sim, tau, true_abundance, true_detection_label) {

  data.frame(
    simulation = sim,
    true_detection_function = true_detection_label,
    true_abundance = true_abundance,

    est_half = NA,
    est_haz = NA,
    est_exp = NA,

    diff_true_half = NA,
    diff_true_haz = NA,
    diff_true_exp = NA,

    diff_half_haz = NA,
    diff_half_exp = NA,
    diff_haz_exp = NA,

    rel_bias_half = NA,
    rel_bias_haz = NA,
    rel_bias_exp = NA,

    availability = availability,
    perceptability = NA,
    tau = tau,
    b = b,
    phi = phi,
    stringsAsFactors = FALSE
  )

}

cat("Starting simulations (parallel over tau)\n")

# function to run all simulations for a single tau value
run_tau_block <- function(tau_value) {
  cat("\n=== Running simulations for tau =", tau_value, "===\n")

  # start with empty results having the right columns
  results_tau <- results[0, ]

  # bSims distance/detection objects use its internal distance units (consistent
  # with rint = rbr and the truncation radius = max(rbr)).
  # Therefore, tau_value should be passed directly in those same internal units.
  tau_m <- tau_value

  for (sim in 1:n_simulations) {

    if (sim %% 50 == 0) cat("tau", tau_value, "- sim", sim, "... ")

    l <- bsims_init(extent = extent, road = 0, edge = 0)
    p <- bsims_populate(l, density = density)
    e <- bsims_animate(p, vocal_rate = phi, duration = duration)

    true_abundance <- get_abundance(p)

    siteCovs <- data.frame(site = 1:n_sites)

    # HALF-NORMAL TRUE
    tryCatch({

      x_list <- vector("list", n_sites)
      d_list <- vector("list", n_sites)

      for (s in 1:n_sites) {

        d <- bsims_detect(e, tau = tau_m,
                          event_type = "both",
                          direction = TRUE)

        d_list[[s]] <- d

        x_list[[s]] <- bsims_transcribe(
          d,
          tint = tbr,
          rint = rbr,
          condition = "det1"
        )

      }

      y <- prepare_for_unmarked(x_list)

      perceptability_value <- calculate_perceptability(
        d_list,
        true_abundance
      )

      results_tau <- rbind(
        results_tau,
        fit_all_keyfuns(
          y,
          true_abundance,
          sim,
          tau_value,
          "half-normal",
          siteCovs,
          perceptability_value
        )
      )

    }, error = function(err) {
      results_tau <<- rbind(results_tau,
                            na_result_row(sim, tau_value, true_abundance, "half-normal"))
    })


    # HAZARD TRUE
    tryCatch({

      x_list <- vector("list", n_sites)
      d_list <- vector("list", n_sites)

      for (s in 1:n_sites) {

        d <- bsims_detect(
          e,
          tau = tau_m,
          dist_fun = function(d, tau, b = 1) {
            1 - exp(-(d / tau)^(-b))
          },
          event_type = "both",
          direction = TRUE
        )

        d_list[[s]] <- d

        x_list[[s]] <- bsims_transcribe(
          d,
          tint = tbr,
          rint = rbr,
          condition = "det1"
        )

      }

      y <- prepare_for_unmarked(x_list)

      perceptability_value <- calculate_perceptability(
        d_list,
        true_abundance
      )

      results_tau <- rbind(
        results_tau,
        fit_all_keyfuns(
          y,
          true_abundance,
          sim,
          tau_value,
          "hazard-rate",
          siteCovs,
          perceptability_value
        )
      )

    }, error = function(err) {
      results_tau <<- rbind(results_tau,
                            na_result_row(sim, tau_value, true_abundance, "hazard-rate"))
    })


    # NEGATIVE EXPONENTIAL TRUE
    tryCatch({

      x_list <- vector("list", n_sites)
      d_list <- vector("list", n_sites)

      for (s in 1:n_sites) {

        d <- bsims_detect(
          e,
          tau = tau_m,
          dist_fun = function(d, tau) exp(-(d / tau)),
          event_type = "both",
          direction = TRUE
        )

        d_list[[s]] <- d

        x_list[[s]] <- bsims_transcribe(
          d,
          tint = tbr,
          rint = rbr,
          condition = "det1"
        )

      }

      y <- prepare_for_unmarked(x_list)

      perceptability_value <- calculate_perceptability(
        d_list,
        true_abundance
      )

      results_tau <- rbind(
        results_tau,
        fit_all_keyfuns(
          y,
          true_abundance,
          sim,
          tau_value,
          "negative-exponential",
          siteCovs,
          perceptability_value
        )
      )

    }, error = function(err) {
      results_tau <<- rbind(results_tau,
                            na_result_row(sim, tau_value, true_abundance,
                                          "negative-exponential"))
    })

  }

  cat("\nCompleted", n_simulations, "simulations for tau =", tau_value, "\n")
  results_tau
}

# run one tau per core (up to length(tau_values))
n_cores <- min(length(tau_values), detectCores())
cat("Using", n_cores, "cores\n")

results_list <- mclapply(tau_values, run_tau_block, mc.cores = n_cores)

results <- do.call(rbind, results_list)

cat("\n=== All simulations complete ===\n")

write.csv(
  results,
  "/Volumes/T7/Chapter2/DataIntegration/Results/DataSimulation/PointCounts25/density_simulation_results_25counts.csv",
  row.names = FALSE
)
