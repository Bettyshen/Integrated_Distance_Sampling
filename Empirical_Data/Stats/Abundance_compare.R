#---
#title: "Compare abundance of species with half-normal and species-specific key function IDS models (with fixed truncation distance 500 m) + 
# Calculate the percentage difference between half-normal and species-specific model"
#author: "Betty Shen"
#date: "2026-01-14" 
#output: html_document
#---
#Goal: the goal of the script is to compare abundance of species with half-normal and species-specific key function IDS models (with fixed truncation distance 500 m)
rm(list = ls())
.libPaths(c("/nfs/stak/users/shenf/.conda/envs/r_env/lib/R/library"))
#Load library
library(unmarked)
library(dplyr)
library(raster)
library(sf)
library(terra)

# ========= Retrieve species list (full + restricted species) ========= #
# Key function
key <- read.csv("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Data/Best_key_functions_trunc_cropped.csv")
# All species list
All_species <- unique(key$Common_Name)
print(All_species)

# Load capping values
cap <- read.csv("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Data/bird_cap_counts.csv")
head(cap)
# Define species-specific capping values for visualization
# Create a named vector with species names as keys and cap values as values
# Only include species that have valid cap values (not NA)
species_caps <- setNames(cap$Count, cap$Common_Name)
# Remove NA values - these species should not be processed
species_caps <- species_caps[!is.na(species_caps)]
cat("Species with cap values:", length(species_caps), "\n")
print(names(species_caps))

# Create a list of species to skip (those with NA cap values)
species_to_skip <- cap$Common_Name[is.na(cap$Count)]
cat("Species to skip (NA cap values):", length(species_to_skip), "\n")
print(species_to_skip) 

# Load species inclusion list
species_include <- read.csv("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Data/species_include.csv")
head(species_include)
# Get species that should be included (Include == "Yes")
species_to_include <- species_include$Species[species_include$Include == "Yes"]
cat("Species to include in study:", length(species_to_include), "\n")
print(species_to_include)


# Function to convert species names between formats
# Convert "Acorn Woodpecker" to "Acorn_Woodpecker" for file names
species_to_filename <- function(species_name) {
  gsub(" ", "_", species_name)
}

# Convert "Acorn_Woodpecker" back to "Acorn Woodpecker" for display
filename_to_species <- function(filename) {
  gsub("_", " ", filename)
}

# =============================================================== #
# ============= Paths to prediction data ============= #
# =============================================================== #

# Paths to prediction directories
    # Half-normal key function IDS models
fixdist_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/PredictionMaps/ContinuVars/FixDistance/Halfnorm"
    # Species-specific key function IDS models
fixdist_species_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/PredictionMaps/ContinuVars/FixDistance/SpeciesSpecif"

# Output directory for abundance comparison
output_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Results/AbundanceComparison"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# =============================================================== #
# ============= Function to extract abundance values ============= #
# =============================================================== #

# Function to extract abundance values from prediction data
extract_abundance <- function(preds) {
  # Check which prediction column exists
  if ("density_pred" %in% names(preds)) {
    return(preds$density_pred)
  } else if ("Predicted" %in% names(preds)) {
    return(preds$Predicted)
  } else {
    # Find first numeric column that's not x or y
    numeric_cols <- names(preds)[sapply(preds, is.numeric)]
    pred_col <- numeric_cols[!numeric_cols %in% c("x", "y")][1]
    if (!is.na(pred_col)) {
      return(preds[[pred_col]])
    } else {
      stop("Could not find prediction column in data")
    }
  }
}

# =============================================================== #
# ============= Function to load predictions for a species ============= #
# =============================================================== #

# Function to load predictions for a species (handles both restricted and non-restricted)
load_predictions <- function(sp, pred_dir) {
  sp_filename <- species_to_filename(sp)
  
  # Try restricted file first (for restricted species)
  restricted_file <- file.path(pred_dir, paste0("RestrictedRaster_1km_Preds_", sp_filename, ".rds"))
  full_file <- file.path(pred_dir, paste0("FullRaster_1km_Preds_", sp_filename, ".rds"))
  
  if (file.exists(restricted_file)) {
    preds <- readRDS(restricted_file)
    return(preds)
  } else if (file.exists(full_file)) {
    preds <- readRDS(full_file)
    return(preds)
  } else {
    return(NULL)
  }
}

# =============================================================== #
# ============= Get list of species to process ============= #
# =============================================================== #

# Get all species that have predictions in at least one directory
all_species_to_process <- c()

# Check FixDistance directory (Half-normal)
fixdist_files <- list.files(fixdist_dir, pattern = "\\.rds$")
fixdist_species <- unique(gsub("(FullRaster_1km_Preds_|RestrictedRaster_1km_Preds_)", "", 
                                gsub("\\.rds$", "", fixdist_files)))
fixdist_species <- filename_to_species(fixdist_species)

# Check FixDistance directory (Species-specific)
fixdist_species_files <- list.files(fixdist_species_dir, pattern = "\\.rds$")
fixdist_specific_species <- unique(gsub("(FullRaster_1km_Preds_|RestrictedRaster_1km_Preds_)", "", 
                                       gsub("\\.rds$", "", fixdist_species_files)))
fixdist_specific_species <- filename_to_species(fixdist_specific_species)

# Get species that have predictions in all four directories
all_species_to_process <- Reduce(intersect, list(fixdist_species, fixdist_specific_species))

# Remove species that should be skipped
all_species_to_process <- all_species_to_process[!all_species_to_process %in% species_to_skip]

# Filter to only include species marked as "Yes" in the inclusion list
all_species_to_process <- all_species_to_process[all_species_to_process %in% species_to_include]

cat("Found", length(all_species_to_process), "species with predictions in all four directories and marked for inclusion:\n")
print(all_species_to_process)

# Sort species for consistent ordering
all_species_to_process <- sort(all_species_to_process)

# =============================================================== #
# ============= Calculate abundance statistics ============= #
# =============================================================== #

# Function to calculate statistics from prediction data frame
# Now extracts LCI and UCI directly from RDS files if available
calculate_stats <- function(preds) {
  # Extract abundance values
  abundance_values <- extract_abundance(preds)
  abundance_values <- abundance_values[!is.na(abundance_values)]
  
  if (length(abundance_values) == 0 || all(abundance_values == 0)) {
    return(list(
      LCI = 0,
      mean = 0,
      UCI = 0,
      median = 0,
      sum = 0
    ))
  }
  
  # Calculate basic statistics
  mean_val <- mean(abundance_values)
  median_val <- median(abundance_values)
  sum_val <- sum(abundance_values)
  
  # Extract LCI and UCI directly from RDS file if available
  # Check for lower/upper columns (from updated prediction files)
  if ("lower" %in% names(preds) && "upper" %in% names(preds)) {
    # Calculate median of lower and upper confidence intervals across all predictions
    # This gives the median LCI and UCI values for the species
    LCI <- median(preds$lower, na.rm = TRUE)
    UCI <- median(preds$upper, na.rm = TRUE)
  } else {
    # If lower/upper columns don't exist, set to NA
    cat("  WARNING: lower/upper columns not found, setting LCI/UCI to NA\n")
    LCI <- NA
    UCI <- NA
  }
  
  return(list(
    LCI = LCI,
    mean = mean_val,
    UCI = UCI,
    median = median_val,
    sum = sum_val
  ))
}

# Initialize results data frame
results_list <- list()

# Process each species
cat("\nProcessing species and calculating abundance statistics...\n")
for (i in seq_along(all_species_to_process)) {
  sp <- all_species_to_process[i]
  cat("\nProcessing species", i, "of", length(all_species_to_process), ":", sp, "\n")
  
  # Load predictions for half-normal model
  preds_halfnorm <- load_predictions(sp, fixdist_dir)
  if (is.null(preds_halfnorm)) {
    cat("  WARNING: Could not load half-normal predictions for", sp, "\n")
    next
  }
  
  # Load predictions for species-specific model
  preds_species <- load_predictions(sp, fixdist_species_dir)
  if (is.null(preds_species)) {
    cat("  WARNING: Could not load species-specific predictions for", sp, "\n")
    next
  }
  
  # Calculate statistics for half-normal model (pass full prediction data frame)
  stats_halfnorm <- calculate_stats(preds_halfnorm)
  
  # Calculate statistics for species-specific model (pass full prediction data frame)
  stats_species <- calculate_stats(preds_species)
  
  # Create row for results
  result_row <- data.frame(
    Species = sp,
    # Half-normal model statistics
    HalfNorm_LCI = stats_halfnorm$LCI,
    HalfNorm_Mean = stats_halfnorm$mean,
    HalfNorm_UCI = stats_halfnorm$UCI,
    HalfNorm_Median = stats_halfnorm$median,
    HalfNorm_Sum = stats_halfnorm$sum,
    # Species-specific model statistics
    SpeciesSpec_LCI = stats_species$LCI,
    SpeciesSpec_Mean = stats_species$mean,
    SpeciesSpec_UCI = stats_species$UCI,
    SpeciesSpec_Median = stats_species$median,
    SpeciesSpec_Sum = stats_species$sum,
    # Percentage difference: (SpeciesSpec_Sum - HalfNorm_Sum) / HalfNorm_Sum
    Percentage_Difference = ifelse(stats_halfnorm$sum != 0, 
                                   ((stats_species$sum - stats_halfnorm$sum) / stats_halfnorm$sum) * 100, 
                                   NA),
    stringsAsFactors = FALSE
  )
  
  results_list[[i]] <- result_row
  
  cat("  Half-normal - Mean:", round(stats_halfnorm$mean, 4), 
      "Sum:", round(stats_halfnorm$sum, 2), "\n")
  cat("  Species-specific - Mean:", round(stats_species$mean, 4), 
      "Sum:", round(stats_species$sum, 2), "\n")
}

# Combine all results
if (length(results_list) > 0) {
  results_df <- do.call(rbind, results_list)
  
  # Export to CSV
  output_file <- file.path(output_dir, "Abundance_Statistics_Comparison.csv")
  write.csv(results_df, output_file, row.names = FALSE)
  cat("\n\nResults exported to:", output_file, "\n")
  cat("Total species processed:", nrow(results_df), "\n")
  
  # Print summary
  cat("\nSummary of results:\n")
  print(head(results_df))
} else {
  cat("\nWARNING: No results to export. Check that prediction files exist for the species.\n")
}
