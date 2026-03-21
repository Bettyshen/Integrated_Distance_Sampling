#Date: 2025-11-17
  
#Goal: Calculate abundance estimate differences by using half-normal vs. species-specific detection key function via 
#data integration (eBird + Oregon 2020) from zero-filled data for multiple species
# Truncation distance is fixed at 500 m

#- Fit models with continuous environmental variables
#- Extract coefficient estimates for predictors
#- Multiple clusters to run integrated distance sampling models
#- Export fitted IDS model 

#=== Load library path for conda environment ===#
.libPaths(c("/nfs/stak/users/shenf/.conda/envs/r_env/lib/R/library"))

#=== Load library & Data ====#
library(TMB)
library(unmarked)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(data.table)

# Base directory
base_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Data"

# Oregon 2020 
OR2020 <- as.data.frame(fread(sprintf("%s/Oregon2020_BirdEnvs_zf_Trunc_cropped_Nov5.csv", base_dir)))

# eBird
eBird <- as.data.frame(fread(sprintf("%s/eBird_zf_cropped_countCap_Nov5.csv", base_dir)))

# Key function
key <- read.csv(sprintf("%s/Best_key_functions_trunc_cropped.csv", base_dir))

#======== Rename predictor variables ========#
# Examine predictor variables
names(OR2020)
# Rename predictor variables
OR2020 <- OR2020 %>%
rename(latitude = Latitude, 
  longitude = Longitude
)


#======== Scale environmental covariates ==========#
# CRITICAL FIX: Calculate scaling parameters from FULL OREGON RASTER
# This ensures consistent scaling between training and prediction
# The full Oregon raster represents all pixels where predictions will be made
cat("Loading full Oregon raster to calculate scaling parameters...\n")
full_raster <- readRDS("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Raster/xy2017_EnvironmentalRaster_Nov5.rds")


names(full_raster)
# List predictors we want to scale
names(OR2020)
names(eBird)
predictors <- c("EVI_median","SR_B1_median","SR_B5_median","ppt_median", "tmean_median", "elevation", "slope", "aspect")

# Select predictor columns from full Oregon raster
# These should match the columns in OR2020 and eBird
common_predictors <- intersect(intersect(predictors, names(full_raster)), 
                               intersect(names(OR2020), names(eBird)))
cat("Calculating scaling parameters from FULL OREGON RASTER...\n")
cat("Predictors to scale:", paste(common_predictors, collapse=", "), "\n")
cat("Full raster has", nrow(full_raster), "pixels\n")

# Extract predictor columns from full Oregon raster
full_raster_predictors <- full_raster[, common_predictors, drop=FALSE]

# Calculate mean and SD for each predictor from FULL OREGON RASTER
scaling_params <- data.frame(
  predictor = common_predictors,
  mean = sapply(full_raster_predictors, function(x) mean(x, na.rm=TRUE)),
  sd = sapply(full_raster_predictors, function(x) sd(x, na.rm=TRUE)),
  stringsAsFactors = FALSE
)

# Print scaling parameters
cat("Scaling parameters (mean, sd) from FULL OREGON RASTER:\n")
print(scaling_params)

# Save scaling parameters for use in prediction script
scaling_params_file <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Data/scaling_parameters_ContinVars.csv"
write.csv(scaling_params, scaling_params_file, row.names = FALSE)
cat("Saved scaling parameters to:", scaling_params_file, "\n")
cat("These parameters will be used for both training data (OR2020, eBird) and prediction data.\n")

# Function to apply consistent scaling using the full raster parameters
apply_consistent_scaling <- function(data, scaling_params) {
  
  for (i in 1:nrow(scaling_params)) {
    pred_name <- scaling_params$predictor[i]
    scaled_name <- paste0(pred_name, "_scaled")
    if (pred_name %in% names(data)) {
      mean_val <- scaling_params$mean[i]
      sd_val <- scaling_params$sd[i]
      
      # Check for zero or very small SD
      if (is.na(sd_val) || sd_val == 0 || abs(sd_val) < 1e-10) {
        cat("WARNING: SD for", pred_name, "is", sd_val, "- using center-only scaling\n")
        data[[scaled_name]] <- data[[pred_name]] - mean_val
      } else {
        data[[scaled_name]] <- (data[[pred_name]] - mean_val) / sd_val
      }
    }
  }
  return(data)
}

# Apply consistent scaling to both training datasets using the FULL RASTER parameters
cat("Applying scaling (from full Oregon raster) to OR2020...\n")
oregon2020_scaled <- apply_consistent_scaling(OR2020, scaling_params)
names(oregon2020_scaled)

cat("Applying scaling (from full Oregon raster) to eBird...\n")
eBird_scaled <- apply_consistent_scaling(eBird, scaling_params)
names(eBird_scaled)

# Clean up full raster from memory (no longer needed)
rm(full_raster, full_raster_predictors)
gc()
#======== Balance dataset in advance =================#
# - Detections vs. non-detections
# OPTIMIZED VERSION for large datasets
balance_bird_occurrence <- function(dataset, detect_col) {
  set.seed(123)
  
  # Convert to data.table for faster operations
  if (!inherits(dataset, "data.table")) {
    dataset <- data.table::as.data.table(dataset)
  }
  
  # Get unique species
  species <- unique(dataset$Common_Name)
  balanced_list <- list()
  
  cat("Processing", length(species), "species for balancing...\n")
  
  for (i in seq_along(species)) {
    sp <- species[i]
    
    # Progress indicator
    if (i %% 10 == 0) {
      cat("Processing species", i, "of", length(species), ":", sp, "\n")
    }
    
    # Use data.table syntax for faster filtering
    sp_data <- dataset[Common_Name == sp]
    
    detect_data <- sp_data[get(detect_col) > 0]
    non_detect_data <- sp_data[get(detect_col) == 0]
    
    n_detect <- nrow(detect_data)
    n_non_detect <- nrow(non_detect_data)
    
    # Only balance if there are both detections and non-detections
    if (n_detect > 0 && n_non_detect > 0) {
      n_sample <- n_detect
      
      if (n_non_detect > n_sample) {
        # Sample without replacement
        sample_idx <- sample(n_non_detect, n_sample)
        non_detect_sample <- non_detect_data[sample_idx]
      } else {
        # Use all non-detections if we have fewer than needed
        non_detect_sample <- non_detect_data
      }
      
      # Combine detections and sampled non-detections
      sp_balanced <- rbind(detect_data, non_detect_sample)
      balanced_list[[sp]] <- sp_balanced
    }
  }
  
  cat("Balancing complete. Combining results...\n")
  
  # Combine all balanced species data
  if (length(balanced_list) > 0) {
    balanced_data <- data.table::rbindlist(balanced_list, fill = TRUE)
    return(as.data.frame(balanced_data))
  } else {
    return(data.frame())
  }
}

# Balance Oregon 2020 & eBird
cat("Starting Oregon 2020 balancing...\n")
OR2020_balanced <- balance_bird_occurrence(oregon2020_scaled, "Occur")
cat("Oregon 2020 balancing complete. Rows:", nrow(OR2020_balanced), "\n")

cat("Starting eBird balancing...\n")
cat("eBird dataset size:", nrow(eBird_scaled), "rows\n")
eBird_balanced <- balance_bird_occurrence(eBird_scaled, "observation_count")
cat("eBird balancing complete. Rows:", nrow(eBird_balanced), "\n")

#======== Compare abundance estimates via species-specific key function & half-normal ======#
#Loop through each species and export coefficient estimates using parallel 

# All species (unique)
unique_SP <- unique(OR2020_balanced$Common_Name)

# Scaled predictor variables
predictors_scaled <- c("EVI_median_scaled","SR_B1_median_scaled",
                       "SR_B5_median_scaled","ppt_median_scaled", 
                       "tmean_median_scaled", "elevation_scaled",
                       "slope_scaled", "aspect_scaled")

# Helper function for failed models
make_na_coef_df <- function(sp, predictors, det_fun){
  data.frame(
    Common_Name = sp,
    Predictor = predictors,
    Coefficient = NA_real_,
    DetectionFunc = det_fun,
    stringsAsFactors = FALSE
  )
}

# Create progress tracking file
progress_file <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Results/progress_tracking_Fixed.txt"
cat("Starting model fitting at", as.character(Sys.time()), "\n", file = progress_file, append = FALSE)
cat("Total species to process:", length(unique_SP), "\n", file = progress_file, append = TRUE)

# Initialize results storage
all_abundance_results <- list()
all_coef_results <- list()

# Setup parallel backend
num_cores <- 5
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterExport(cl, varlist = c("OR2020_balanced", "eBird_balanced", "key", "make_na_coef_df", "predictors_scaled"))

clusterEvalQ(cl, {
  library(unmarked)
  library(dplyr)
  library(tidyr)
})

# ========== MODIFIED PARALLEL LOOP WITH IMMEDIATE SAVING ========== #
cat("Starting parallel processing of", length(unique_SP), "species...\n")

# Process species in batches to allow for immediate saving
batch_size <- 10  # Process 10 species at a time
total_species <- length(unique_SP)

for (batch_start in seq(1, total_species, by = batch_size)) {
  batch_end <- min(batch_start + batch_size - 1, total_species)
  batch_species <- unique_SP[batch_start:batch_end]
  
  cat("Processing batch", ceiling(batch_start/batch_size), "of", ceiling(total_species/batch_size), 
      "- Species", batch_start, "to", batch_end, "\n")
  
  # Parallel loop for current batch
  batch_results <- foreach(sp = batch_species, .packages = c("unmarked", "dplyr", "tidyr")) %dopar% {
    set.seed(123)
    
    cat("Processing species:", sp, "\n")
    
    ##====== Oregon 2020 data =====##
    sp_data <- OR2020_balanced %>% filter(Common_Name == sp)
    unique_sp_data <- sp_data %>%
      group_by(Unique_Observation_ID) %>%
      slice_sample(n = 1) %>%
      ungroup() %>% as.data.frame()
    
    trunc_dist <- key %>% filter(Common_Name == sp) %>% pull(Truncation_Distance)
    if (length(trunc_dist) == 0) return(NULL)
    
    unique_sp_data$Unique_Observation_ID <- as.factor(unique_sp_data$Unique_Observation_ID)
    dist_breaks <- seq(0, trunc_dist + 10, by = 15)
    dis_sp <- formatDistData(unique_sp_data, distCol = "Distance",
                             transectNameCol = "Unique_Observation_ID",
                             dist.breaks = dist_breaks)
    
    covs_sp <- as.data.frame(unique_sp_data[, c("Unique_Observation_ID", predictors_scaled)])
    umf_sp <- unmarkedFrameDS(y = as.matrix(dis_sp), siteCovs = covs_sp, 
                              survey = "point", dist.breaks = dist_breaks, unitsIn = "m")
    
    ##====== eBird data ======##
    sp_eB <- eBird_balanced %>% filter(Common_Name == sp)
    sp_eB$sampling_event_identifier <- as.factor(sp_eB$sampling_event_identifier)
    covs_eB <- as.data.frame(sp_eB[, c("sampling_event_identifier", predictors_scaled)])
    Mpc <- nrow(sp_eB)
    db2 <- c(0, 500)
    umf_eB <- unmarkedFrameDS(y = matrix(sp_eB$observation_count, Mpc, 1), 
                              siteCovs = covs_eB, survey = "point", dist.breaks = db2, unitsIn = "m")
    umf_pc <- unmarkedFramePCount(y = umf_eB@y, siteCovs = umf_eB@siteCovs)
    
    # === Best key function of the species === #
    key_function_raw <- key %>% filter(Common_Name == sp) %>% pull(Best_Key_Function)
    
    # Map key function names to unmarked format
    key_function <- case_when(
      key_function_raw == "half_normal" ~ "halfnorm",
      key_function_raw == "negative_exponential" ~ "exp", 
      key_function_raw == "hazard_rate" ~ "hazard",
      TRUE ~ key_function_raw  # fallback to original if no match
    )
    
    # Fit models
    IDS.model <- tryCatch({
      IDS(lambdaformula = ~ EVI_median_scaled + SR_B1_median_scaled + SR_B5_median_scaled +
            ppt_median_scaled + tmean_median_scaled + elevation_scaled + slope_scaled + aspect_scaled,
          dataDS = umf_sp, dataPC = umf_pc, keyfun = key_function,
          unitsOut = "kmsq", maxDistPC = 500)
    }, error = function(e) {
      cat("Error fitting species-specific model for", sp, ":", e$message, "\n")
      NULL
    })
    
    IDS.half <- tryCatch({
      IDS(lambdaformula = ~ EVI_median_scaled + SR_B1_median_scaled + SR_B5_median_scaled +
            ppt_median_scaled + tmean_median_scaled + elevation_scaled + slope_scaled + aspect_scaled,
          dataDS = umf_sp, dataPC = umf_pc, keyfun = "halfnorm",
          unitsOut = "kmsq", maxDistPC = 500)
    }, error = function(e) {
      cat("Error fitting half-normal model for", sp, ":", e$message, "\n")
      NULL
    })
    
    # Prepare species name for file saving
    sp_clean <- gsub("/", "", sp)
    
    # ========== SAVE MODELS DIRECTLY IN PARALLEL WORKER ========== #
    # Save models immediately if they exist (in parallel worker to avoid serialization issues)
    if (!is.null(IDS.model)) {
      tryCatch({
        file_path_spec <- paste0("/nfs/stak/users/shenf/hpc-share/IntegrationModel/TrainedModel/ContinuVars/FixDistance/SpeciesSpecif/",
                                sp_clean, ".rds")
        saveRDS(IDS.model, file = file_path_spec)
        # Force file system sync
        flush.connection(stdout())
        cat("Saved species-specific model for:", sp, "\n")
      }, error = function(e) {
        cat("ERROR saving species-specific model for", sp, ":", e$message, "\n")
      })
    }
    
    if (!is.null(IDS.half)) {
      tryCatch({
        file_path_half <- paste0("/nfs/stak/users/shenf/hpc-share/IntegrationModel/TrainedModel/ContinuVars/FixDistance/Halfnorm/",
                                sp_clean, ".rds")
        saveRDS(IDS.half, file = file_path_half)
        # Force file system sync
        flush.connection(stdout())
        cat("Saved half-normal model for:", sp, "\n")
      }, error = function(e) {
        cat("ERROR saving half-normal model for", sp, ":", e$message, "\n")
      })
    }
    
    # Add NA if the model was not trained successfully
    if (is.null(IDS.model) || is.null(IDS.half)) {
      abundance_out <- data.frame(Common_Name = sp,
                                  Species_Specific_Abundance = NA,
                                  Species_Specific_LCIabundance = NA,
                                  Species_Specific_UCIabundance = NA,
                                  Half_Normal_Abundance = NA,
                                  Half_Normal_LCIabundance = NA,
                                  Half_Normal_UCIabundance = NA,
                                  AIC_Specific = NA,
                                  AIC_HalfNormal = NA)
      
      coefs_out <- do.call(rbind, list(
        if (is.null(IDS.model)) make_na_coef_df(sp, predictors, key_function_raw),
        if (is.null(IDS.half))  make_na_coef_df(sp, predictors, "halfnorm")
      ))
      
      return(list(abundance = abundance_out, coefs = coefs_out))
    }
    
    # Predict abundance on training data itself
    predict.IDS <- predict(IDS.model, "lam", backTransform = TRUE)
    predict.half.IDS <- predict(IDS.half, "lam", backTransform = TRUE)
    
    # Sum all predicted abundance per site
    # Species-specific key function
    # Total abundance
    species_specific_abundance <- sum(predict.IDS$Predicted)
    # Lower bound abundance
    species_specific_LCIabundance <- sum(predict.IDS$lower)
    # Upper bound abundance
    species_specific_UCIabundance <- sum(predict.IDS$upper)
    
    # Half-normal
    # Total abundance
    half_normal_abundance <- sum(predict.half.IDS$Predicted)
    # Lower bound abundance
    half_normal_LCIabundance <- sum(predict.half.IDS$lower)
    # Upper bound abundance
    half_normal_UCIabundance <- sum(predict.half.IDS$upper)
    
    # Extract abundance information both species-specific vs. half-normal key function
    abundance_out <- data.frame(Common_Name = sp,
                                Species_Specific_Abundance = species_specific_abundance,
                                Species_Specific_LCIabundance = species_specific_LCIabundance,
                                Species_Specific_UCIabundance = species_specific_UCIabundance,
                                Half_Normal_Abundance = half_normal_abundance,
                                Half_Normal_LCIabundance = half_normal_LCIabundance,
                                Half_Normal_UCIabundance = half_normal_UCIabundance,
                                AIC_Specific = IDS.model@AIC,
                                AIC_HalfNormal = IDS.half@AIC)
    
    # Extract coefficient estimates for each environmental predictor
    coef_Specfunc <- IDS.model@estimates@estimates[["lam"]]@estimates
    coef_halfnorm <- IDS.half@estimates@estimates[["lam"]]@estimates
    
    df_Specfunc <- data.frame(
      Common_Name = sp,
      Predictor = names(coef_Specfunc),
      Coefficient = unname(coef_Specfunc),
      DetectionFunc = key_function_raw
    )
    
    df_halfnorm <- data.frame(
      Common_Name = sp,
      Predictor = names(coef_halfnorm),
      Coefficient = unname(coef_halfnorm),
      DetectionFunc = "halfnorm"
    )
    
    coefs_out <- rbind(df_Specfunc, df_halfnorm)
    
    # Models already saved above, no need to return them (avoids serialization issues)
    return(list(abundance = abundance_out, coefs = coefs_out))
  }
  
  # ========== IMMEDIATE RESULTS SAVING ========== #
  # Combine results from current batch
  valid_batch_results <- Filter(Negate(is.null), batch_results)
  
  cat("Batch", ceiling(batch_start/batch_size), "- Valid results:", length(valid_batch_results), "out of", length(batch_results), "\n")
  
  # Check how many models were actually saved
  if (length(valid_batch_results) > 0) {
    # Count saved files for this batch
    saved_specific <- sum(sapply(batch_species, function(sp) {
      sp_clean <- gsub("/", "", sp)
      file.exists(paste0("/nfs/stak/users/shenf/hpc-share/IntegrationModel/TrainedModel/ContinuVars/FixDistance/SpeciesSpecif/", sp_clean, ".rds"))
    }))
    saved_half <- sum(sapply(batch_species, function(sp) {
      sp_clean <- gsub("/", "", sp)
      file.exists(paste0("/nfs/stak/users/shenf/hpc-share/IntegrationModel/TrainedModel/ContinuVars/FixDistance/Halfnorm/", sp_clean, ".rds"))
    }))
    cat("Batch", ceiling(batch_start/batch_size), "summary: Found", saved_specific, 
        "species-specific models and", saved_half, "half-normal models saved\n")
  }
  
  if (length(valid_batch_results) > 0) {
    batch_abundance <- do.call(rbind, lapply(valid_batch_results, function(x) x$abundance))
    batch_coefs <- do.call(rbind, lapply(valid_batch_results, function(x) x$coefs))
    
    # Append to overall results
    all_abundance_results[[length(all_abundance_results) + 1]] <- batch_abundance
    all_coef_results[[length(all_coef_results) + 1]] <- batch_coefs
    
    # Save intermediate results
    current_abundance <- do.call(rbind, all_abundance_results)
    current_coefs <- do.call(rbind, all_coef_results)
    
    write.csv(current_abundance, "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Results/Abundance_Comparison_continuVars_FixDistance.csv", row.names = FALSE)
    write.csv(current_coefs, "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Results/IDS_coefficients_AllSp_continuVars_FixDistance.csv", row.names = FALSE)
    
    # Update progress
    cat("Completed batch", ceiling(batch_start/batch_size), "- Species", batch_start, "to", batch_end, 
        "at", as.character(Sys.time()), "\n", file = progress_file, append = TRUE)
    cat("Models saved so far:", nrow(current_abundance), "species\n", file = progress_file, append = TRUE)
  }
}

# Stop parallel backend
stopCluster(cl)

# Final results
final_abundance <- do.call(rbind, all_abundance_results)
final_coefs <- do.call(rbind, all_coef_results)

# Save final results
write.csv(final_abundance, "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Results/Abundance_Comparison_continuVars_FixDistance_FINAL.csv", row.names = FALSE)
write.csv(final_coefs, "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Results/IDS_coefficients_AllSp_continuVars_FixDistance_FINAL.csv", row.names = FALSE)

# View in console
print(final_abundance)
print(final_coefs)

cat("All processing complete at", as.character(Sys.time()), "\n", file = progress_file, append = TRUE)
cat("Total species processed:", nrow(final_abundance), "\n", file = progress_file, append = TRUE)
