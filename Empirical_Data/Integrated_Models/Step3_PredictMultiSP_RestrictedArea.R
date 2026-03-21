#---
#title: "Multi-species Prediction Maps for half normal + species-specific key function IDS model (Restricted + non-restricted) with continuous variables"
#author: "Betty Shen"
#date: "2025-12-08"
#output: html_document
#---
#Goal: Use consistent scaling parameters from full training dataset (fix distance 500 m as truncation distance)
#This ensures consistent scaling across all species and fixes inverted habitat preferences
# Focus on continuous variables: EVI, SR_B1, SR_B5, ppt, tmean, elevation, slope, aspect
# This script predicts abundance using both half-normal and species-specific detection function models

.libPaths(c("/nfs/stak/users/shenf/.conda/envs/r_env/lib/R/library"))
#Load library
library(unmarked)
library(dplyr)
library(raster)
library(sf)
library(terra)

# Load species list
key <- read.csv("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Data/Best_key_functions_trunc_cropped.csv")
All_species <- unique(key$Common_Name)
print(All_species)

# Filter species list to only include those with successfully trained models
model_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/TrainedModel/ContinuVars/FixDistance/Halfnorm"
available_models <- list.files(model_dir, pattern = "\\.rds$")
available_species <- gsub("\\.rds$", "", available_models)
length(available_species)
# Filter All_species to only include those with available models
All_species <- All_species[All_species %in% available_species]



cat("Species with available models (", length(All_species), "total):\n")
print(All_species)

# Function to convert species names between formats
species_to_filename <- function(species_name) {
  gsub(" ", "_", species_name)
}

filename_to_species <- function(filename) {
  gsub("_", " ", filename)
}

# =============================================================== #
# ============= STEP 1: LOAD CONSISTENT SCALING PARAMETERS ============= #
# =============================================================== #

# CRITICAL FIX: Load scaling parameters from training script
# These parameters were calculated from the FULL OREGON RASTER
# The training script calculates mean/SD from all Oregon raster pixels,
# then applies those same parameters to OR2020 and eBird training data
# This ensures predictions use the same scaling as the models were trained on
scaling_params_file <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Data/scaling_parameters_ContinVars.csv"

if (!file.exists(scaling_params_file)) {
  cat("ERROR: Scaling parameters file not found!\n")
  cat("File expected at:", scaling_params_file, "\n")
  cat("Please run the model training script first to generate scaling parameters.\n")
  cat("The training script will save scaling parameters to this file.\n")
  stop("Cannot proceed without scaling parameters file.")
}

cat("Loading scaling parameters from training data...\n")
scaling_params <- read.csv(scaling_params_file, stringsAsFactors = FALSE)
cat("Loaded scaling parameters for", nrow(scaling_params), "predictors:\n")
print(scaling_params)

# Create a lookup for quick access
scaling_lookup <- setNames(scaling_params$mean, scaling_params$predictor)
scaling_sd_lookup <- setNames(scaling_params$sd, scaling_params$predictor)

# Function to apply consistent scaling using the full Oregon rasterdata parameters
apply_consistent_scaling <- function(data) {
  # Apply scaling for each predictor
  for (i in 1:nrow(scaling_params)) {
    pred_name <- scaling_params$predictor[i]
    scaled_name <- paste0(pred_name, "_scaled")
    
    if (pred_name %in% names(data)) {
      mean_val <- scaling_params$mean[i]
      sd_val <- scaling_params$sd[i]
      
      # Check for zero or very small SD (would cause division issues)
      if (is.na(sd_val) || sd_val == 0 || abs(sd_val) < 1e-10) {
        cat("WARNING: SD for", pred_name, "is", sd_val, "- using unscaled values\n")
        data[[scaled_name]] <- data[[pred_name]] - mean_val
      } else {
        data[[scaled_name]] <- (data[[pred_name]] - mean_val) / sd_val
      }
    } else {
      cat("WARNING: Predictor", pred_name, "not found in data. Available columns:", paste(names(data), collapse=", "), "\n")
    }
  }
  return(data)
}

cat("Scaling function ready. Using parameters from the full Oregon raster.\n")

# =============================================================== #
# ============= STEP 2: PREDICT FOR RESTRICTED SPECIES ============= #
# =============================================================== #

# Load Restricted area raster directory
raster_restricted_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/RestrictedArea/RestrictedPredArea/ContinuVars"

# Get restricted species list
raster_restrictedTIFF_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/RestrictedArea/TiffFile"
species_files <- list.files(raster_restrictedTIFF_dir, pattern = "_PredArea.tif$")
restricted_species <- gsub("_PredArea.tif", "", species_files)
cat("Found", length(restricted_species), "restricted species to process:\n")
print(restricted_species)

## We'll read and prepare each species' restricted raster inside the loop below
chunk_size <- 50000

# Path to trained models (half normal detection models)
model_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/TrainedModel/ContinuVars/FixDistance/Halfnorm"
model_paths <- file.path(model_dir, paste0(filename_to_species(restricted_species), ".rds"))

# Read all models into a named list
IDS_models_restricted <- setNames(lapply(model_paths, readRDS), restricted_species)

print("Predicting for restricted species with CONSISTENT scaling...")

for (sp in restricted_species) {
  cat("Predicting for", sp, "...\n")
  IDS.model <- IDS_models_restricted[[sp]]
  # Load raster for this species
  sp_filename <- species_to_filename(sp)
  sp_path <- file.path(raster_restricted_dir, paste0("xy2017_EnvironRaster_", sp_filename, "_FullRaster_Restricted.rds"))
  
  # Check if file exists
  if (!file.exists(sp_path)) {
    cat("Warning: File not found:", sp_path, "\n")
    cat("Skipping species:", sp, "\n")
    next
  }
  
  restricted_raster <- readRDS(sp_path)
  print(names(restricted_raster))

  # Only keep the columns that match the model's covariates
  # Actual structure: columns 3:10 = SR_B1_median, SR_B5_median, EVI_median, ppt_median, 
  #                   tmean_median, elevation, slope, aspect
  newdata_restricted <- restricted_raster[,c(3:10)]
  cat("Column names in newdata_restricted:", paste(names(newdata_restricted), collapse=", "), "\n")
  
  # Check ranges of all predictors before scaling
  cat("\nPre-scaling value ranges:\n")
  for(col in names(newdata_restricted)) {
    if(is.numeric(newdata_restricted[[col]])) {
      col_range <- range(newdata_restricted[[col]], na.rm=TRUE)
      cat(sprintf("  %s: [%.4f, %.4f]\n", col, col_range[1], col_range[2]))
    }
  }

  # Check how many pixels have actual environmental data vs. 0 vs. NAs
  cat("Total pixels in raster:", nrow(newdata_restricted), "\n")

  # Count pixels with actual environmental data vs. non-NA vs. non-zeros
  # Check for non-NA values (not checking sum > 0 since scaled values can be negative/zero)
  has_data <- rowSums(!is.na(newdata_restricted[,1:8])) > 0
  has_zeros <- rowSums(!is.na(newdata_restricted[,1:8])) > 0 & rowSums(newdata_restricted[,1:8], na.rm=TRUE) == 0
  has_nas <- rowSums(!is.na(newdata_restricted[,1:8])) == 0

  cat("Pixels with actual environmental data (non-zero):", sum(has_data), "\n")
  cat("Pixels outside presence area (all zeros):", sum(has_zeros), "\n")
  cat("Pixels with NAs (outside raster extent):", sum(has_nas), "\n")

  # ============= CRITICAL FIX: CONSISTENT SCALING ============= #
  # Apply consistent scaling using the FULL training dataset parameters
  cat("Applying CONSISTENT scaling using full training dataset parameters...\n")
  
  # Create the scaled dataframe with columns that match the model's covariates
  newdata_restricted_scaled <- data.frame(
    newdata_restricted, # original columns
    EVI_median_scaled = rep(0, nrow(newdata_restricted)),
    SR_B1_median_scaled = rep(0, nrow(newdata_restricted)),
    SR_B5_median_scaled = rep(0, nrow(newdata_restricted)),
    ppt_median_scaled = rep(0, nrow(newdata_restricted)),
    tmean_median_scaled = rep(0, nrow(newdata_restricted)),
    elevation_scaled = rep(0, nrow(newdata_restricted)),
    slope_scaled = rep(0, nrow(newdata_restricted)),
    aspect_scaled = rep(0, nrow(newdata_restricted))
  )
  
  # Apply consistent scaling to ALL pixels using the training data parameters
  newdata_restricted_scaled <- apply_consistent_scaling(newdata_restricted_scaled)

  # Use the scaled data for prediction
  newdata_restricted_final <- newdata_restricted_scaled[, c("EVI_median_scaled", "SR_B1_median_scaled", "SR_B5_median_scaled", 
                                                            "ppt_median_scaled", "tmean_median_scaled", "elevation_scaled", 
                                                            "slope_scaled", "aspect_scaled")]
  names(newdata_restricted_final)
  
  # CRITICAL: Validate scaled values are within reasonable ranges
  # Extreme scaled values (>5 or <-5) may indicate scaling mismatch or data issues
  cat("\nValidating scaled values for", sp, ":\n")
  cat("Expected scaled value range: typically -3 to +3 (within 3 SD of mean)\n")
  cat("Values outside -5 to +5 may indicate coordinate system mismatch or data issues\n\n")
  
  for (col in names(newdata_restricted_final)) {
    col_data <- newdata_restricted_final[[col]]
    
    # Get the original predictor name
    orig_col <- gsub("_scaled", "", col)
    
    # Get scaling parameters for this predictor
    mean_val <- NA
    sd_val <- NA
    if(orig_col %in% scaling_params$predictor) {
      param_row <- scaling_params[scaling_params$predictor == orig_col, ]
      mean_val <- param_row$mean
      sd_val <- param_row$sd
      cat(sprintf("%s:\n", col))
      cat(sprintf("  Scaling params: mean=%.4f, sd=%.4f\n", mean_val, sd_val))
    } else {
      cat(sprintf("%s:\n", col))
      cat("  WARNING: No scaling parameters found for this predictor\n")
    }
    
    extreme_high <- sum(col_data > 5, na.rm=TRUE)
    extreme_low <- sum(col_data < -5, na.rm=TRUE)
    col_range <- range(col_data, na.rm=TRUE)
    
    cat(sprintf("  Scaled range: [%.4f, %.4f]\n", col_range[1], col_range[2]))
    
    if (extreme_high > 0 || extreme_low > 0) {
      cat(sprintf("  WARNING: %d values > 5, %d values < -5\n", extreme_high, extreme_low))
      if(!is.na(mean_val) && !is.na(sd_val) && sd_val > 0) {
        if(extreme_high > 0) {
          cat(sprintf("    Max value: %.4f (%.1f SD above mean)\n", col_range[2], (col_range[2] - mean_val) / sd_val))
        }
        if(extreme_low > 0) {
          cat(sprintf("    Min value: %.4f (%.1f SD below mean)\n", col_range[1], (mean_val - col_range[1]) / sd_val))
        }
      }
    }
    
    # Check for all zeros or all NAs
    if (all(col_data == 0, na.rm=TRUE)) {
      cat("  WARNING: All values are 0!\n")
    }
    if (all(is.na(col_data))) {
      cat("  ERROR: All values are NA!\n")
    }
    cat("\n")
  }

  # Prepare chunking based on this species' raster size
  n_rows <- nrow(newdata_restricted_final)
  n_chunks <- ceiling(n_rows / chunk_size)

  all_preds <- list()
  # Predict density per pixel for each species using chunked prediction
  for (i in 1:n_chunks) {
    start_row <- ((i-1) * chunk_size) + 1
    end_row <- min(i * chunk_size, n_rows)
    chunk_data <- newdata_restricted_final[start_row:end_row, ]

    # CRITICAL FIX: Check for pixels with valid environmental data
    # After scaling, row sums can be negative, zero, or positive, so we only check for non-NA values
    # This ensures we have actual environmental data (not all NAs)
    chunk_has_data <- rowSums(!is.na(chunk_data)) > 0

    if (sum(chunk_has_data) > 0) {
      # Predict only for pixels with data
      chunk_data_valid <- chunk_data[chunk_has_data, ]
      chunk_preds <- predict(IDS.model, newdata = chunk_data_valid, type = "lam")

      # CRITICAL FIX: Keep full data frame with Predicted, lower, upper, SE columns
      # unmarked predict() returns a data frame with "Predicted", "lower", "upper", "SE" columns
      if (is.null(chunk_preds)) {
        cat("WARNING: predict() returned NULL for chunk", i, "- creating empty data frame\n")
        chunk_preds <- data.frame(Predicted = rep(0, sum(chunk_has_data)),
                                  SE = rep(0, sum(chunk_has_data)),
                                  lower = rep(0, sum(chunk_has_data)),
                                  upper = rep(0, sum(chunk_has_data)))
      } else if (is.data.frame(chunk_preds)) {
        # Keep the full data frame - it should have Predicted, lower, upper, SE
        # Check if all expected columns exist, if not create them with zeros
        if(!"Predicted" %in% names(chunk_preds)) {
          if("lam" %in% names(chunk_preds)) {
            chunk_preds$Predicted <- chunk_preds$lam
          } else {
            numeric_cols <- names(chunk_preds)[sapply(chunk_preds, is.numeric)]
            if(length(numeric_cols) > 0) {
              chunk_preds$Predicted <- chunk_preds[[numeric_cols[1]]]
              cat("Warning: Using column '", numeric_cols[1], "' as Predicted\n")
            } else {
              cat("ERROR: No numeric column found in predict() output for chunk", i, "\n")
              chunk_preds <- data.frame(Predicted = rep(0, sum(chunk_has_data)),
                                        SE = rep(0, sum(chunk_has_data)),
                                        lower = rep(0, sum(chunk_has_data)),
                                        upper = rep(0, sum(chunk_has_data)))
            }
          }
        }
        # Ensure lower and upper columns exist
        if(!"lower" %in% names(chunk_preds)) {
          chunk_preds$lower <- rep(0, nrow(chunk_preds))
          cat("Warning: 'lower' column not found, setting to 0\n")
        }
        if(!"upper" %in% names(chunk_preds)) {
          chunk_preds$upper <- rep(0, nrow(chunk_preds))
          cat("Warning: 'upper' column not found, setting to 0\n")
        }
        if(!"SE" %in% names(chunk_preds)) {
          chunk_preds$SE <- rep(0, nrow(chunk_preds))
          cat("Warning: 'SE' column not found, setting to 0\n")
        }
      } else if (is.list(chunk_preds)) {
        # Convert list to data frame
        chunk_preds <- as.data.frame(chunk_preds)
        # Ensure required columns exist
        if(!"Predicted" %in% names(chunk_preds)) {
          chunk_preds$Predicted <- unlist(chunk_preds)
        }
        if(!"lower" %in% names(chunk_preds)) chunk_preds$lower <- rep(0, nrow(chunk_preds))
        if(!"upper" %in% names(chunk_preds)) chunk_preds$upper <- rep(0, nrow(chunk_preds))
        if(!"SE" %in% names(chunk_preds)) chunk_preds$SE <- rep(0, nrow(chunk_preds))
      } else {
        cat("WARNING: Unexpected predict() output format for chunk", i, "- converting\n")
        chunk_preds <- data.frame(Predicted = as.numeric(chunk_preds),
                                  SE = rep(0, length(chunk_preds)),
                                  lower = rep(0, length(chunk_preds)),
                                  upper = rep(0, length(chunk_preds)))
      }
      
      # Validate that chunk_preds is a data frame with required columns
      if(!is.data.frame(chunk_preds)) {
        cat("ERROR: chunk_preds is not a data frame after processing for chunk", i, "\n")
        chunk_preds <- data.frame(Predicted = rep(0, sum(chunk_has_data)),
                                  SE = rep(0, sum(chunk_has_data)),
                                  lower = rep(0, sum(chunk_has_data)),
                                  upper = rep(0, sum(chunk_has_data)))
      }
      
      # Validate number of rows
      if(nrow(chunk_preds) != sum(chunk_has_data)) {
        cat("WARNING: Row mismatch in chunk", i, "- expected", sum(chunk_has_data), "got", nrow(chunk_preds), "\n")
        if(nrow(chunk_preds) < sum(chunk_has_data)) {
          # Pad with zeros
          n_missing <- sum(chunk_has_data) - nrow(chunk_preds)
          padding <- data.frame(Predicted = rep(0, n_missing),
                                SE = rep(0, n_missing),
                                lower = rep(0, n_missing),
                                upper = rep(0, n_missing))
          chunk_preds <- rbind(chunk_preds, padding)
        } else {
          chunk_preds <- chunk_preds[1:sum(chunk_has_data), ]
        }
      }
      
      # Check for all-zero predictions
      if(all(chunk_preds$Predicted == 0, na.rm=TRUE)) {
        cat("WARNING: All predictions are 0 for chunk", i, "of species", sp, "\n")
        cat("This may indicate a problem with the model or input data.\n")
      }
      
      # Create full prediction data frame with zeros for pixels without data
      full_chunk_preds <- data.frame(
        Predicted = rep(0, nrow(chunk_data)),
        SE = rep(0, nrow(chunk_data)),
        lower = rep(0, nrow(chunk_data)),
        upper = rep(0, nrow(chunk_data))
      )
      full_chunk_preds[chunk_has_data, ] <- chunk_preds

      all_preds[[i]] <- full_chunk_preds
    } else {
      # All pixels in this chunk are outside presence area (all zeros)
      all_preds[[i]] <- data.frame(
        Predicted = rep(0, nrow(chunk_data)),
        SE = rep(0, nrow(chunk_data)),
        lower = rep(0, nrow(chunk_data)),
        upper = rep(0, nrow(chunk_data))
      )
    }

    cat("Completed chunk", i, "of", n_chunks, "- chunk size:", nrow(all_preds[[i]]), "\n")

  }
  # Combine all chunks into a single data frame
  preds <- do.call(rbind, all_preds)
  
  # Debug: Check if preds has the correct number of rows
  cat("Total predictions generated:", nrow(preds), "\n")
  cat("Expected rows (restricted raster rows):", nrow(restricted_raster), "\n")
  
  # CRITICAL: Validate prediction values before saving
  cat("Summary of all predictions for", sp, ":\n")
  cat("Min Predicted:", min(preds$Predicted, na.rm=TRUE), "\n")
  cat("Max Predicted:", max(preds$Predicted, na.rm=TRUE), "\n")
  cat("Mean Predicted:", mean(preds$Predicted, na.rm=TRUE), "\n")
  cat("Number of non-zero predictions:", sum(preds$Predicted > 0, na.rm=TRUE), "\n")
  cat("Number of zero predictions:", sum(preds$Predicted == 0, na.rm=TRUE), "\n")
  cat("Number of NA predictions:", sum(is.na(preds$Predicted)), "\n")
  
  # Check if all predictions are zero
  if(all(preds$Predicted == 0, na.rm=TRUE) || (sum(preds$Predicted > 0, na.rm=TRUE) == 0)) {
    cat("ERROR: All prediction values are 0 for", sp, "!\n")
    cat("This indicates a serious problem. Possible causes:\n")
    cat("  1. The predict() function is not returning valid predictions\n")
    cat("  2. The model is not working correctly\n")
    cat("  3. The input data format is incorrect\n")
    cat("  4. The condition for 'chunk_has_data' is too strict\n")
    cat("Please check the model and input data for this species.\n")
  }
  
  # Ensure preds has the correct number of rows
  if (nrow(preds) != nrow(restricted_raster)) {
    cat("ERROR: Prediction data frame row mismatch!\n")
    cat("This will cause issues in the visualization script.\n")
    # Try to fix by padding with zeros or truncating
    if (nrow(preds) < nrow(restricted_raster)) {
      cat("Padding predictions with zeros...\n")
      n_missing <- nrow(restricted_raster) - nrow(preds)
      padding <- data.frame(Predicted = rep(0, n_missing),
                            SE = rep(0, n_missing),
                            lower = rep(0, n_missing),
                            upper = rep(0, n_missing))
      preds <- rbind(preds, padding)
    } else {
      cat("Truncating predictions to match raster length...\n")
      preds <- preds[1:nrow(restricted_raster), ]
    }
  }

  # Create final prediction dataframe with coordinates and all prediction columns
  final_preds <- data.frame(
    x = restricted_raster$x,
    y = restricted_raster$y,
    density_pred = preds$Predicted,
    SE = preds$SE,
    lower = preds$lower,
    upper = preds$upper
  )
  
  # Final validation
  cat("Final prediction dataframe summary:\n")
  cat("Number of rows:", nrow(final_preds), "\n")
  cat("Number of non-zero density_pred:", sum(final_preds$density_pred > 0, na.rm=TRUE), "\n")
  cat("density_pred range:", range(final_preds$density_pred, na.rm=TRUE), "\n")
  cat("Mean density:", mean(final_preds$density_pred, na.rm=TRUE), "\n")
  cat("Median density:", median(final_preds$density_pred, na.rm=TRUE), "\n")
  cat("Columns saved:", paste(names(final_preds), collapse=", "), "\n")
  
  # Save predictions
  pred_file <- file.path("/nfs/stak/users/shenf/hpc-share/IntegrationModel/PredictionMaps/ContinuVars/FixDistance/Halfnorm", 
                         paste0("RestrictedRaster_1km_Preds_", gsub(' ', '_', sp), ".rds"))
  saveRDS(final_preds, pred_file)
  cat("Done predicting for", sp, "...\n")
}


# =============================================================== #
# ============= STEP 3: PREDICT FOR NON-RESTRICTED SPECIES ============= #
# =============================================================== #

# Get species list that are not in restricted area
not_restricted_species <- setdiff(All_species, restricted_species)
cat("Found", length(not_restricted_species), "non-restricted species to process:\n")
print(not_restricted_species)

# Load raster for not restricted species
raster <- readRDS("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Raster/xy2017_EnvironmentalRaster_Nov5.rds")

print(names(raster))

# Only keep the columns that match the model's covariates
# Actual structure: columns 3:10 = SR_B1_median, SR_B5_median, EVI_median, ppt_median, 
#                   tmean_median, elevation, slope, aspect
names(raster)
newdata <- raster[,c(3:10)]

# Check ranges of all predictors before scaling
cat("\nPre-scaling value ranges:\n")
for(col in names(newdata)) {
  if(is.numeric(newdata[[col]])) {
    col_range <- range(newdata[[col]], na.rm=TRUE)
    cat(sprintf("  %s: [%.4f, %.4f]\n", col, col_range[1], col_range[2]))
  }
}




# Apply consistent scaling using the same parameters
cat("\nApplying consistent scaling to non-restricted species...\n")
newdata_scaled <- apply_consistent_scaling(newdata)

# Use the scaled data for prediction
newdata <- newdata_scaled[, c("EVI_median_scaled", "SR_B1_median_scaled", "SR_B5_median_scaled", 
                             "ppt_median_scaled", "tmean_median_scaled", "elevation_scaled", 
                             "slope_scaled", "aspect_scaled")]
names(newdata)
summary(newdata)
# CRITICAL: Validate scaled values are within reasonable ranges
# Extreme scaled values (>5 or <-5) may indicate scaling mismatch or data issues
cat("\nValidating scaled values for non-restricted species:\n")
cat("Expected scaled value range: typically -3 to +3 (within 3 SD of mean)\n")
cat("Values outside -5 to +5 may indicate coordinate system mismatch or data issues\n\n")

for (col in names(newdata)) {
  col_data <- newdata[[col]]
  
  # Get the original predictor name
  orig_col <- gsub("_scaled", "", col)
  
  # Get scaling parameters for this predictor
  mean_val <- NA
  sd_val <- NA
  if(orig_col %in% scaling_params$predictor) {
    param_row <- scaling_params[scaling_params$predictor == orig_col, ]
    mean_val <- param_row$mean
    sd_val <- param_row$sd
    cat(sprintf("%s:\n", col))
    cat(sprintf("  Scaling params: mean=%.4f, sd=%.4f\n", mean_val, sd_val))
  } else {
    cat(sprintf("%s:\n", col))
    cat("  WARNING: No scaling parameters found for this predictor\n")
  }
  
  extreme_high <- sum(col_data > 5, na.rm=TRUE)
  extreme_low <- sum(col_data < -5, na.rm=TRUE)
  col_range <- range(col_data, na.rm=TRUE)
  
  cat(sprintf("  Scaled range: [%.4f, %.4f]\n", col_range[1], col_range[2]))
  
  if (extreme_high > 0 || extreme_low > 0) {
    cat(sprintf("  WARNING: %d values > 5, %d values < -5\n", extreme_high, extreme_low))
    if(!is.na(mean_val) && !is.na(sd_val) && sd_val > 0) {
      if(extreme_high > 0) {
        cat(sprintf("    Max value: %.4f (%.1f SD above mean)\n", col_range[2], (col_range[2] - mean_val) / sd_val))
      }
      if(extreme_low > 0) {
        cat(sprintf("    Min value: %.4f (%.1f SD below mean)\n", col_range[1], (mean_val - col_range[1]) / sd_val))
      }
    }
  }
  
  # Check for all zeros or all NAs
  if (all(col_data == 0, na.rm=TRUE)) {
    cat("  WARNING: All values are 0!\n")
  }
  if (all(is.na(col_data))) {
    cat("  ERROR: All values are NA!\n")
  }
  cat("\n")
}

print("Load Multi-species half normal IDS model")

# Path to trained models (half normal detection models)
model_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/TrainedModel/ContinuVars/FixDistance/Halfnorm"
model_paths <- file.path(model_dir, paste0(not_restricted_species, ".rds"))

# Read all models into a named list
IDS_models <- setNames(lapply(model_paths, readRDS), not_restricted_species)

# Predict density per pixel for each species using chunked prediction
chunk_size <- 50000
n_rows <- nrow(newdata)
n_chunks <- ceiling(n_rows / chunk_size)

# Predict density per pixel for selected species
print("Predicting for non-restricted species with CONSISTENT scaling...")

for (sp in not_restricted_species) {
  cat("Predicting for", sp, "...\n")
  IDS.model <- IDS_models[[sp]]
  all_preds <- list()
  for (i in 1:n_chunks) {
    start_row <- ((i-1) * chunk_size) + 1
    end_row <- min(i * chunk_size, n_rows)
    chunk_data <- newdata[start_row:end_row, ]
    chunk_preds <- predict(IDS.model, newdata = chunk_data, type = "lam")
    
    # CRITICAL FIX: Handle different output formats from predict()
    if (is.null(chunk_preds)) {
      cat("WARNING: predict() returned NULL for chunk", i, "- creating empty data frame\n")
      chunk_preds <- data.frame(Predicted = rep(0, nrow(chunk_data)))
    } else if (!is.data.frame(chunk_preds)) {
      # Convert to data frame if it's a vector or matrix
      if(is.vector(chunk_preds) || is.numeric(chunk_preds)) {
        chunk_preds <- data.frame(Predicted = as.numeric(chunk_preds))
      } else if(is.matrix(chunk_preds)) {
        if("Predicted" %in% colnames(chunk_preds)) {
          chunk_preds <- data.frame(Predicted = chunk_preds[, "Predicted"])
        } else {
          chunk_preds <- data.frame(Predicted = chunk_preds[, 1])
        }
      } else {
        cat("WARNING: Unexpected predict() output format for chunk", i, "- converting\n")
        chunk_preds <- data.frame(Predicted = as.numeric(chunk_preds))
      }
    }
    
    # Validate chunk predictions
    if("Predicted" %in% names(chunk_preds)) {
      if(all(chunk_preds$Predicted == 0, na.rm=TRUE)) {
        cat("WARNING: All predictions are 0 for chunk", i, "of species", sp, "\n")
      }
    }
    
    all_preds[[i]] <- chunk_preds
    cat("Completed chunk", i, "of", n_chunks, "\n")
  }
  preds <- do.call(rbind, all_preds)
  
  # CRITICAL: Validate prediction values before saving
  if("Predicted" %in% names(preds)) {
    cat("Summary of all predictions for", sp, ":\n")
    cat("Min:", min(preds$Predicted, na.rm=TRUE), "\n")
    cat("Max:", max(preds$Predicted, na.rm=TRUE), "\n")
    cat("Mean:", mean(preds$Predicted, na.rm=TRUE), "\n")
    cat("Median:", median(preds$Predicted, na.rm=TRUE), "\n")
    cat("Number of non-zero predictions:", sum(preds$Predicted > 0, na.rm=TRUE), "\n")
    cat("Number of zero predictions:", sum(preds$Predicted == 0, na.rm=TRUE), "\n")
    cat("Number of NA predictions:", sum(is.na(preds$Predicted)), "\n")
    
    # Check if all predictions are zero
    if(all(preds$Predicted == 0, na.rm=TRUE) || (sum(preds$Predicted > 0, na.rm=TRUE) == 0)) {
      cat("ERROR: All prediction values are 0 for", sp, "!\n")
      cat("This indicates a serious problem. Please check the model and input data.\n")
    }
  } else {
    cat("WARNING: 'Predicted' column not found in predictions for", sp, "\n")
    cat("Available columns:", paste(names(preds), collapse=", "), "\n")
  }
  
  # Save predictions
  pred_file <- file.path("/nfs/stak/users/shenf/hpc-share/IntegrationModel/PredictionMaps/ContinuVars/FixDistance/Halfnorm", 
                         paste0("FullRaster_1km_Preds_", gsub(' ', '_', sp), ".rds"))
  saveRDS(preds, pred_file)
  cat("Done predicting for", sp, "...\n")
}

cat("All species predictions completed with consistent scaling!\n")

# =============================================================== #
# ============= STEP 4: PREDICT FOR RESTRICTED SPECIES (SPECIES-SPECIFIC) ============= #
# =============================================================== #

# Filter species list to only include those with species-specific models available
model_dir_spec <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/TrainedModel/ContinuVars/FixDistance/SpeciesSpecif"
available_models_spec <- list.files(model_dir_spec, pattern = "\\.rds$")
available_species_spec <- gsub("\\.rds$", "", available_models_spec)
cat("Species with available species-specific models (", length(available_species_spec), "total):\n")
print(available_species_spec)

# Get restricted species that have species-specific models
# Convert restricted_species (filenames) to species names for comparison
restricted_species_names <- filename_to_species(restricted_species)
restricted_species_spec <- restricted_species[restricted_species_names %in% available_species_spec]
cat("Found", length(restricted_species_spec), "restricted species with species-specific models to process:\n")
print(restricted_species_spec)

# Path to trained models (species-specific detection models)
model_paths_spec_restricted <- file.path(model_dir_spec, paste0(filename_to_species(restricted_species_spec), ".rds"))

# Read all models into a named list
IDS_models_restricted_spec <- setNames(lapply(model_paths_spec_restricted, readRDS), restricted_species_spec)

print("Predicting for restricted species with species-specific models and CONSISTENT scaling...")

for (sp in restricted_species_spec) {
  cat("Predicting for", sp, "using species-specific model...\n")
  IDS.model <- IDS_models_restricted_spec[[sp]]
  # Load raster for this species
  sp_filename <- species_to_filename(sp)
  sp_path <- file.path(raster_restricted_dir, paste0("xy2017_EnvironRaster_", sp_filename, "_FullRaster_Restricted.rds"))
  
  # Check if file exists
  if (!file.exists(sp_path)) {
    cat("Warning: File not found:", sp_path, "\n")
    cat("Skipping species:", sp, "\n")
    next
  }
  
  restricted_raster <- readRDS(sp_path)
  print(names(restricted_raster))

  # Only keep the columns that match the model's covariates
  newdata_restricted <- restricted_raster[,c(3:10)]
  cat("Column names in newdata_restricted:", paste(names(newdata_restricted), collapse=", "), "\n")
  
  # Check ranges of all predictors before scaling
  cat("\nPre-scaling value ranges:\n")
  for(col in names(newdata_restricted)) {
    if(is.numeric(newdata_restricted[[col]])) {
      col_range <- range(newdata_restricted[[col]], na.rm=TRUE)
      cat(sprintf("  %s: [%.4f, %.4f]\n", col, col_range[1], col_range[2]))
    }
  }

  # Check how many pixels have actual environmental data vs. 0 vs. NAs
  cat("Total pixels in raster:", nrow(newdata_restricted), "\n")
  has_data <- rowSums(!is.na(newdata_restricted[,1:8])) > 0
  has_zeros <- rowSums(!is.na(newdata_restricted[,1:8])) > 0 & rowSums(newdata_restricted[,1:8], na.rm=TRUE) == 0
  has_nas <- rowSums(!is.na(newdata_restricted[,1:8])) == 0

  cat("Pixels with actual environmental data (non-zero):", sum(has_data), "\n")
  cat("Pixels outside presence area (all zeros):", sum(has_zeros), "\n")
  cat("Pixels with NAs (outside raster extent):", sum(has_nas), "\n")

  # ============= CRITICAL FIX: CONSISTENT SCALING ============= #
  cat("Applying CONSISTENT scaling using full training dataset parameters...\n")
  
  # Create the scaled dataframe with columns that match the model's covariates
  newdata_restricted_scaled <- data.frame(
    newdata_restricted, # original columns
    EVI_median_scaled = rep(0, nrow(newdata_restricted)),
    SR_B1_median_scaled = rep(0, nrow(newdata_restricted)),
    SR_B5_median_scaled = rep(0, nrow(newdata_restricted)),
    ppt_median_scaled = rep(0, nrow(newdata_restricted)),
    tmean_median_scaled = rep(0, nrow(newdata_restricted)),
    elevation_scaled = rep(0, nrow(newdata_restricted)),
    slope_scaled = rep(0, nrow(newdata_restricted)),
    aspect_scaled = rep(0, nrow(newdata_restricted))
  )
  
  # Apply consistent scaling to ALL pixels using the training data parameters
  newdata_restricted_scaled <- apply_consistent_scaling(newdata_restricted_scaled)

  # Use the scaled data for prediction
  newdata_restricted_final <- newdata_restricted_scaled[, c("EVI_median_scaled", "SR_B1_median_scaled", "SR_B5_median_scaled", 
                                                            "ppt_median_scaled", "tmean_median_scaled", "elevation_scaled", 
                                                            "slope_scaled", "aspect_scaled")]
  names(newdata_restricted_final)
  
  # CRITICAL: Validate scaled values are within reasonable ranges
  cat("\nValidating scaled values for", sp, "(species-specific model):\n")
  cat("Expected scaled value range: typically -3 to +3 (within 3 SD of mean)\n")
  cat("Values outside -5 to +5 may indicate coordinate system mismatch or data issues\n\n")
  
  for (col in names(newdata_restricted_final)) {
    col_data <- newdata_restricted_final[[col]]
    orig_col <- gsub("_scaled", "", col)
    mean_val <- NA
    sd_val <- NA
    if(orig_col %in% scaling_params$predictor) {
      param_row <- scaling_params[scaling_params$predictor == orig_col, ]
      mean_val <- param_row$mean
      sd_val <- param_row$sd
      cat(sprintf("%s:\n", col))
      cat(sprintf("  Scaling params: mean=%.4f, sd=%.4f\n", mean_val, sd_val))
    } else {
      cat(sprintf("%s:\n", col))
      cat("  WARNING: No scaling parameters found for this predictor\n")
    }
    
    extreme_high <- sum(col_data > 5, na.rm=TRUE)
    extreme_low <- sum(col_data < -5, na.rm=TRUE)
    col_range <- range(col_data, na.rm=TRUE)
    
    cat(sprintf("  Scaled range: [%.4f, %.4f]\n", col_range[1], col_range[2]))
    
    if (extreme_high > 0 || extreme_low > 0) {
      cat(sprintf("  WARNING: %d values > 5, %d values < -5\n", extreme_high, extreme_low))
      if(!is.na(mean_val) && !is.na(sd_val) && sd_val > 0) {
        if(extreme_high > 0) {
          cat(sprintf("    Max value: %.4f (%.1f SD above mean)\n", col_range[2], (col_range[2] - mean_val) / sd_val))
        }
        if(extreme_low > 0) {
          cat(sprintf("    Min value: %.4f (%.1f SD below mean)\n", col_range[1], (mean_val - col_range[1]) / sd_val))
        }
      }
    }
    
    if (all(col_data == 0, na.rm=TRUE)) {
      cat("  WARNING: All values are 0!\n")
    }
    if (all(is.na(col_data))) {
      cat("  ERROR: All values are NA!\n")
    }
    cat("\n")
  }

  # Prepare chunking based on this species' raster size
  n_rows <- nrow(newdata_restricted_final)
  n_chunks <- ceiling(n_rows / chunk_size)

  all_preds <- list()
  # Predict density per pixel for each species using chunked prediction
  for (i in 1:n_chunks) {
    start_row <- ((i-1) * chunk_size) + 1
    end_row <- min(i * chunk_size, n_rows)
    chunk_data <- newdata_restricted_final[start_row:end_row, ]

    # CRITICAL FIX: Check for pixels with valid environmental data
    chunk_has_data <- rowSums(!is.na(chunk_data)) > 0

    if (sum(chunk_has_data) > 0) {
      # Predict only for pixels with data
      chunk_data_valid <- chunk_data[chunk_has_data, ]
      chunk_preds <- predict(IDS.model, newdata = chunk_data_valid, type = "lam")

      # CRITICAL FIX: Keep full data frame with Predicted, lower, upper, SE columns
      # unmarked predict() returns a data frame with "Predicted", "lower", "upper", "SE" columns
      if (is.null(chunk_preds)) {
        cat("WARNING: predict() returned NULL for chunk", i, "- creating empty data frame\n")
        chunk_preds <- data.frame(Predicted = rep(0, sum(chunk_has_data)),
                                  SE = rep(0, sum(chunk_has_data)),
                                  lower = rep(0, sum(chunk_has_data)),
                                  upper = rep(0, sum(chunk_has_data)))
      } else if (is.data.frame(chunk_preds)) {
        # Keep the full data frame - it should have Predicted, lower, upper, SE
        # Check if all expected columns exist, if not create them with zeros
        if(!"Predicted" %in% names(chunk_preds)) {
          if("lam" %in% names(chunk_preds)) {
            chunk_preds$Predicted <- chunk_preds$lam
          } else {
            numeric_cols <- names(chunk_preds)[sapply(chunk_preds, is.numeric)]
            if(length(numeric_cols) > 0) {
              chunk_preds$Predicted <- chunk_preds[[numeric_cols[1]]]
              cat("Warning: Using column '", numeric_cols[1], "' as Predicted\n")
            } else {
              cat("ERROR: No numeric column found in predict() output for chunk", i, "\n")
              chunk_preds <- data.frame(Predicted = rep(0, sum(chunk_has_data)),
                                        SE = rep(0, sum(chunk_has_data)),
                                        lower = rep(0, sum(chunk_has_data)),
                                        upper = rep(0, sum(chunk_has_data)))
            }
          }
        }
        # Ensure lower and upper columns exist
        if(!"lower" %in% names(chunk_preds)) {
          chunk_preds$lower <- rep(0, nrow(chunk_preds))
          cat("Warning: 'lower' column not found, setting to 0\n")
        }
        if(!"upper" %in% names(chunk_preds)) {
          chunk_preds$upper <- rep(0, nrow(chunk_preds))
          cat("Warning: 'upper' column not found, setting to 0\n")
        }
        if(!"SE" %in% names(chunk_preds)) {
          chunk_preds$SE <- rep(0, nrow(chunk_preds))
          cat("Warning: 'SE' column not found, setting to 0\n")
        }
      } else if (is.list(chunk_preds)) {
        # Convert list to data frame
        chunk_preds <- as.data.frame(chunk_preds)
        # Ensure required columns exist
        if(!"Predicted" %in% names(chunk_preds)) {
          chunk_preds$Predicted <- unlist(chunk_preds)
        }
        if(!"lower" %in% names(chunk_preds)) chunk_preds$lower <- rep(0, nrow(chunk_preds))
        if(!"upper" %in% names(chunk_preds)) chunk_preds$upper <- rep(0, nrow(chunk_preds))
        if(!"SE" %in% names(chunk_preds)) chunk_preds$SE <- rep(0, nrow(chunk_preds))
      } else {
        cat("WARNING: Unexpected predict() output format for chunk", i, "- converting\n")
        chunk_preds <- data.frame(Predicted = as.numeric(chunk_preds),
                                  SE = rep(0, length(chunk_preds)),
                                  lower = rep(0, length(chunk_preds)),
                                  upper = rep(0, length(chunk_preds)))
      }
      
      # Validate that chunk_preds is a data frame with required columns
      if(!is.data.frame(chunk_preds)) {
        cat("ERROR: chunk_preds is not a data frame after processing for chunk", i, "\n")
        chunk_preds <- data.frame(Predicted = rep(0, sum(chunk_has_data)),
                                  SE = rep(0, sum(chunk_has_data)),
                                  lower = rep(0, sum(chunk_has_data)),
                                  upper = rep(0, sum(chunk_has_data)))
      }
      
      # Validate number of rows
      if(nrow(chunk_preds) != sum(chunk_has_data)) {
        cat("WARNING: Row mismatch in chunk", i, "- expected", sum(chunk_has_data), "got", nrow(chunk_preds), "\n")
        if(nrow(chunk_preds) < sum(chunk_has_data)) {
          # Pad with zeros
          n_missing <- sum(chunk_has_data) - nrow(chunk_preds)
          padding <- data.frame(Predicted = rep(0, n_missing),
                                SE = rep(0, n_missing),
                                lower = rep(0, n_missing),
                                upper = rep(0, n_missing))
          chunk_preds <- rbind(chunk_preds, padding)
        } else {
          chunk_preds <- chunk_preds[1:sum(chunk_has_data), ]
        }
      }
      
      # Check for all-zero predictions
      if(all(chunk_preds$Predicted == 0, na.rm=TRUE)) {
        cat("WARNING: All predictions are 0 for chunk", i, "of species", sp, "\n")
        cat("This may indicate a problem with the model or input data.\n")
      }
      
      # Create full prediction data frame with zeros for pixels without data
      full_chunk_preds <- data.frame(
        Predicted = rep(0, nrow(chunk_data)),
        SE = rep(0, nrow(chunk_data)),
        lower = rep(0, nrow(chunk_data)),
        upper = rep(0, nrow(chunk_data))
      )
      full_chunk_preds[chunk_has_data, ] <- chunk_preds

      all_preds[[i]] <- full_chunk_preds
    } else {
      # All pixels in this chunk are outside presence area (all zeros)
      all_preds[[i]] <- data.frame(
        Predicted = rep(0, nrow(chunk_data)),
        SE = rep(0, nrow(chunk_data)),
        lower = rep(0, nrow(chunk_data)),
        upper = rep(0, nrow(chunk_data))
      )
    }

    cat("Completed chunk", i, "of", n_chunks, "- chunk size:", nrow(all_preds[[i]]), "\n")

  }
  # Combine all chunks into a single data frame
  preds <- do.call(rbind, all_preds)
  
  # Debug: Check if preds has the correct number of rows
  cat("Total predictions generated:", nrow(preds), "\n")
  cat("Expected rows (restricted raster rows):", nrow(restricted_raster), "\n")
  
  # CRITICAL: Validate prediction values before saving
  cat("Summary of all predictions for", sp, "(species-specific model):\n")
  cat("Min Predicted:", min(preds$Predicted, na.rm=TRUE), "\n")
  cat("Max Predicted:", max(preds$Predicted, na.rm=TRUE), "\n")
  cat("Mean Predicted:", mean(preds$Predicted, na.rm=TRUE), "\n")
  cat("Number of non-zero predictions:", sum(preds$Predicted > 0, na.rm=TRUE), "\n")
  cat("Number of zero predictions:", sum(preds$Predicted == 0, na.rm=TRUE), "\n")
  cat("Number of NA predictions:", sum(is.na(preds$Predicted)), "\n")
  
  # Check if all predictions are zero
  if(all(preds$Predicted == 0, na.rm=TRUE) || (sum(preds$Predicted > 0, na.rm=TRUE) == 0)) {
    cat("ERROR: All prediction values are 0 for", sp, "!\n")
    cat("This indicates a serious problem. Possible causes:\n")
    cat("  1. The predict() function is not returning valid predictions\n")
    cat("  2. The model is not working correctly\n")
    cat("  3. The input data format is incorrect\n")
    cat("  4. The condition for 'chunk_has_data' is too strict\n")
    cat("Please check the model and input data for this species.\n")
  }
  
  # Ensure preds has the correct number of rows
  if (nrow(preds) != nrow(restricted_raster)) {
    cat("ERROR: Prediction data frame row mismatch!\n")
    cat("This will cause issues in the visualization script.\n")
    if (nrow(preds) < nrow(restricted_raster)) {
      cat("Padding predictions with zeros...\n")
      n_missing <- nrow(restricted_raster) - nrow(preds)
      padding <- data.frame(Predicted = rep(0, n_missing),
                            SE = rep(0, n_missing),
                            lower = rep(0, n_missing),
                            upper = rep(0, n_missing))
      preds <- rbind(preds, padding)
    } else {
      cat("Truncating predictions to match raster length...\n")
      preds <- preds[1:nrow(restricted_raster), ]
    }
  }

  # Create final prediction dataframe with coordinates and all prediction columns
  final_preds <- data.frame(
    x = restricted_raster$x,
    y = restricted_raster$y,
    density_pred = preds$Predicted,
    SE = preds$SE,
    lower = preds$lower,
    upper = preds$upper
  )
  
  # Final validation
  cat("Final prediction dataframe summary:\n")
  cat("Number of rows:", nrow(final_preds), "\n")
  cat("Number of non-zero density_pred:", sum(final_preds$density_pred > 0, na.rm=TRUE), "\n")
  cat("density_pred range:", range(final_preds$density_pred, na.rm=TRUE), "\n")
  cat("Mean density:", mean(final_preds$density_pred, na.rm=TRUE), "\n")
  cat("Median density:", median(final_preds$density_pred, na.rm=TRUE), "\n")
  cat("Columns saved:", paste(names(final_preds), collapse=", "), "\n")
  
  # Save predictions to SpeciesSpecif directory
  pred_file <- file.path("/nfs/stak/users/shenf/hpc-share/IntegrationModel/PredictionMaps/ContinuVars/FixDistance/SpeciesSpecif", 
                         paste0("RestrictedRaster_1km_Preds_", gsub(' ', '_', sp), ".rds"))
  saveRDS(final_preds, pred_file)
  cat("Done predicting for", sp, "with species-specific model...\n")
}


# =============================================================== #
# ============= STEP 5: PREDICT FOR NON-RESTRICTED SPECIES (SPECIES-SPECIFIC) ============= #
# =============================================================== #

# Get species list that are not in restricted area and have species-specific models
# Convert restricted_species_spec (filenames) to species names for comparison
restricted_species_spec_names <- filename_to_species(restricted_species_spec)
not_restricted_species_spec <- setdiff(available_species_spec, restricted_species_spec_names)
cat("Found", length(not_restricted_species_spec), "non-restricted species with species-specific models to process:\n")
print(not_restricted_species_spec)

# Reuse the same scaled newdata from Step 3 (already scaled and validated)
cat("\nUsing already-scaled data from Step 3 for species-specific predictions...\n")

print("Load Multi-species species-specific IDS model")

# Path to trained models (species-specific detection models)
model_paths_spec <- file.path(model_dir_spec, paste0(not_restricted_species_spec, ".rds"))

# Read all models into a named list
IDS_models_spec <- setNames(lapply(model_paths_spec, readRDS), not_restricted_species_spec)

# Predict density per pixel for each species using chunked prediction
chunk_size <- 50000
n_rows <- nrow(newdata)
n_chunks <- ceiling(n_rows / chunk_size)

# Predict density per pixel for selected species
print("Predicting for non-restricted species with species-specific models and CONSISTENT scaling...")

for (sp in not_restricted_species_spec) {
  cat("Predicting for", sp, "using species-specific model...\n")
  IDS.model <- IDS_models_spec[[sp]]
  all_preds <- list()
  for (i in 1:n_chunks) {
    start_row <- ((i-1) * chunk_size) + 1
    end_row <- min(i * chunk_size, n_rows)
    chunk_data <- newdata[start_row:end_row, ]
    chunk_preds <- predict(IDS.model, newdata = chunk_data, type = "lam")
    
    # CRITICAL FIX: Handle different output formats from predict()
    if (is.null(chunk_preds)) {
      cat("WARNING: predict() returned NULL for chunk", i, "- creating empty data frame\n")
      chunk_preds <- data.frame(Predicted = rep(0, nrow(chunk_data)))
    } else if (!is.data.frame(chunk_preds)) {
      # Convert to data frame if it's a vector or matrix
      if(is.vector(chunk_preds) || is.numeric(chunk_preds)) {
        chunk_preds <- data.frame(Predicted = as.numeric(chunk_preds))
      } else if(is.matrix(chunk_preds)) {
        if("Predicted" %in% colnames(chunk_preds)) {
          chunk_preds <- data.frame(Predicted = chunk_preds[, "Predicted"])
        } else {
          chunk_preds <- data.frame(Predicted = chunk_preds[, 1])
        }
      } else {
        cat("WARNING: Unexpected predict() output format for chunk", i, "- converting\n")
        chunk_preds <- data.frame(Predicted = as.numeric(chunk_preds))
      }
    }
    
    # Validate chunk predictions
    if("Predicted" %in% names(chunk_preds)) {
      if(all(chunk_preds$Predicted == 0, na.rm=TRUE)) {
        cat("WARNING: All predictions are 0 for chunk", i, "of species", sp, "\n")
      }
    }
    
    all_preds[[i]] <- chunk_preds
    cat("Completed chunk", i, "of", n_chunks, "\n")
  }
  preds <- do.call(rbind, all_preds)
  
  # CRITICAL: Validate prediction values before saving
  if("Predicted" %in% names(preds)) {
    cat("Summary of all predictions for", sp, "(species-specific model):\n")
    cat("Min:", min(preds$Predicted, na.rm=TRUE), "\n")
    cat("Max:", max(preds$Predicted, na.rm=TRUE), "\n")
    cat("Mean:", mean(preds$Predicted, na.rm=TRUE), "\n")
    cat("Median:", median(preds$Predicted, na.rm=TRUE), "\n")
    cat("Number of non-zero predictions:", sum(preds$Predicted > 0, na.rm=TRUE), "\n")
    cat("Number of zero predictions:", sum(preds$Predicted == 0, na.rm=TRUE), "\n")
    cat("Number of NA predictions:", sum(is.na(preds$Predicted)), "\n")
    
    # Check if all predictions are zero
    if(all(preds$Predicted == 0, na.rm=TRUE) || (sum(preds$Predicted > 0, na.rm=TRUE) == 0)) {
      cat("ERROR: All prediction values are 0 for", sp, "!\n")
      cat("This indicates a serious problem. Please check the model and input data.\n")
    }
  } else {
    cat("WARNING: 'Predicted' column not found in predictions for", sp, "\n")
    cat("Available columns:", paste(names(preds), collapse=", "), "\n")
  }
  
  # Save predictions to SpeciesSpecif directory
  pred_file <- file.path("/nfs/stak/users/shenf/hpc-share/IntegrationModel/PredictionMaps/ContinuVars/FixDistance/SpeciesSpecif", 
                         paste0("FullRaster_1km_Preds_", gsub(' ', '_', sp), ".rds"))
  saveRDS(preds, pred_file)
  cat("Done predicting for", sp, "with species-specific model...\n")
}

cat("All species predictions completed with consistent scaling for both half-normal and species-specific models!\n")
