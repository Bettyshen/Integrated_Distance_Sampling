# Goal: export pixel values along with x, y coordinates, restricted to species presence areas
# Added new covariates: slope, aspect
# Date: 2025-11-09
rm(list=ls())
#==== For All species ====#
# Add library path for conda environment
.libPaths(c("/nfs/stak/users/shenf/.conda/envs/r_env/lib/R/library"))
# Load files & library
library(dplyr)
library(terra)

# ========================== Load Oregon raster ======================# 
raster <- rast("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Raster/EnvironmentalRaster_Nov5.tif")


# === Reproject raster to Oregon State Plane EPSG: 2992 === #
new_crs <- "EPSG:2992"
cat("Reprojecting to:", new_crs, "\n")
reprojected_raster <- project(raster, new_crs)

# =============== Load restricted area for All species ============== #
raster_restricted_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/RestrictedArea/TiffFile"
# Get species list
species_files <- list.files(raster_restricted_dir, pattern = "_PredArea.tif$")
species <- gsub("_PredArea.tif", "", species_files)
cat("Found", length(species), "species to process:\n")
print(species)

# =============== Export restricted area with xy coordinates ============== #
# =============== Process each species individually ============== #
# Create the export directory if it doesn't exist
export_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/RestrictedArea/RestrictedPredArea/ContinuVars"
if(!dir.exists(export_dir)) {
  dir.create(export_dir, recursive = TRUE)
}

# Process each species one by one
for(i in 1:length(species)) {
  sp <- species[i]
  cat("\n=== Processing species", i, "of", length(species), ":", sp, "===\n")
  
  # Load restricted area for this specific species
  restricted_file <- file.path(raster_restricted_dir, paste0(sp, "_PredArea.tif"))
  
  # Check if file exists
  if(!file.exists(restricted_file)) {
    cat("Warning: File not found:", restricted_file, "\n")
    next
  }
  
  tryCatch({
    # Load the restricted area raster for this species
    restricted <- rast(restricted_file)
    cat("Loaded restricted area for", sp, "\n")
    
    # If no CRS is set, we need to assign one EPSG:2992
    if(is.na(crs(restricted)) || crs(restricted) == "") {
      cat("No CRS detected. Assigning (EPSG:2992) to restricted raster...\n")
      crs(restricted) <- "EPSG:2992"
    }
    
    # Reproject restricted area to match environmental raster CRS (2992(Oregon State Plane))
    restricted_reprojected <- project(restricted, new_crs)
    
    # === Check restricted area CRS === #
    cat("Restricted area CRS:", crs(restricted), "\n")
    cat("Restricted area dimensions:", dim(restricted), "\n")
    
    # === Ensure both rasters have the same CRS and resolution === #
    cat("Ensuring CRS and resolution compatibility...\n")
    
    # === First, reproject the restricted area to match the environmental raster to make sure they have the same CRS=== #
    cat("Reprojecting restricted area to match environmental raster CRS...\n")
    restricted_reprojected <- project(restricted, crs(reprojected_raster), method = "near")
    
    # === Resample to match the environmental raster's resolution and extent === #
    cat("Resampling restricted area to match environmental raster resolution...\n")
    restricted_resampled <- resample(restricted_reprojected, reprojected_raster, method = "near")
    
    cat("CRS and resolution alignment completed.\n")
    
    # ============= Create a mask where 1 = presence area, 0 = outside presence area ============= #
    # === Create a mask where 1 = presence area, 0 = outside presence area === #
    cat("Creating presence area mask...\n")
    presence_mask <- restricted_resampled
    presence_mask[!is.na(presence_mask) & presence_mask > 0] <- 1
    presence_mask[is.na(presence_mask) | presence_mask == 0] <- 0
    
    cat("Mask created. Presence pixels:", sum(presence_mask[], na.rm = TRUE), "\n")
    cat("Outside presence pixels:", sum(presence_mask[] == 0, na.rm = TRUE), "\n")
    
    # === Apply the mask to the environmental raster === #
    # This will set all values outside the presence area to 0
    cat("Applying mask to environmental raster...\n")
    restricted_raster <- reprojected_raster * presence_mask
    cat("Masking completed.\n")
    
    # ============= Convert restricted raster to dataframe with coordinates ============= #
    # === Convert restricted raster to dataframe with coordinates === #
    # Create two versions: one with only presence area data, one with full raster extent
    
    # Version 1: Only cells with actual environmental data (for analysis)
    # === Version 1: Only cells with actual environmental data (for analysis) === #
    cat("Extracting only cells with environmental data...\n")
    data_cells <- which(!is.na(restricted_raster[]) & rowSums(restricted_raster[] != 0, na.rm = TRUE) > 0)
    cat("Found", length(data_cells), "cells with environmental data\n")
    
    if(length(data_cells) > 0) {
      # === Get coordinates for cells with data === #
      coords <- xyFromCell(restricted_raster, data_cells)
      
      # === Get environmental values for these cells === #
      env_values <- restricted_raster[data_cells]
      
      # === Create dataframe with only the data cells === #
      df_restricted <- data.frame(
        x = coords[,1],
        y = coords[,2],
        longitude = coords[,1],
        env_values
      )
      
      cat("Created restricted dataframe with", nrow(df_restricted), "rows with environmental data\n")
      
    } else {
      cat("No cells with environmental data found!\n")
      df_restricted <- data.frame()
    }
    
    # Version 2: Full raster extent with zeros outside presence area (for visualization/prediction)
    # === Version 2: Full raster extent with zeros outside presence area (for visualization/prediction) === #
    cat("Creating full raster version with zeros outside presence area...\n")
    df_full_restricted <- as.data.frame(restricted_raster, xy = TRUE, na.rm = FALSE)
    
    # Add longitude column (same as x coordinate)
    # The prediction script expects columns 3:11 to include: EVI_median, SR_B1_median, SR_B5_median, 
    # ppt_median, tmean_median, elevation, longitude, slope, aspect (9 columns total)
    df_full_restricted$longitude <- df_full_restricted$x
    
    # Print column information for debugging
    cat("Column names in df_full_restricted:", paste(names(df_full_restricted), collapse = ", "), "\n")
    cat("Number of columns:", ncol(df_full_restricted), "\n")
    cat("Columns 3:11 will be:", paste(names(df_full_restricted)[3:min(11, ncol(df_full_restricted))], collapse = ", "), "\n")
    
    # === Check the full raster version === #
    cat("Full raster dataframe dimensions:", dim(df_full_restricted), "\n")
    cat("Full raster non-zero environmental cells:", sum(rowSums(df_full_restricted[,3:ncol(df_full_restricted)] != 0, na.rm = TRUE) > 0), "\n")
    cat("Full raster zero/NA cells:", sum(rowSums(df_full_restricted[,3:ncol(df_full_restricted)] == 0, na.rm = TRUE) == ncol(df_full_restricted)-2), "\n")
    
    # === Create a clean version with only the presence area data (this should be the same as df_restricted now) === #
    df_restricted_clean <- df_restricted
    
    # === Export restricted dataframes === #
    # === Version 1: Only cells with actual environmental in restricted area (for analysis) === #
    restricted_file_out <- file.path(export_dir, 
                               paste0("xy2017_EnvironRaster_", gsub(' ', '_', sp), "_Restricted.rds"))
    # === Version 2: Full raster extent with zeros outside presence area (for visualization/prediction) === #
    full_restricted_file <- file.path(export_dir, 
                               paste0("xy2017_EnvironRaster_", gsub(' ', '_', sp), "_FullRaster_Restricted.rds"))
    
    saveRDS(df_restricted_clean, restricted_file_out)
    saveRDS(df_full_restricted, full_restricted_file)
    
    cat("Successfully processed and saved data for", sp, "\n")
    cat("Files saved:\n")
    cat("  -", restricted_file_out, "\n")
    cat("  -", full_restricted_file, "\n")
    
    # Clean up memory
    rm(restricted, restricted_reprojected, restricted_resampled, presence_mask, 
       restricted_raster, df_restricted, df_full_restricted, df_restricted_clean)
    gc()
    
  }, error = function(e) {
    cat("Error processing", sp, ":", e$message, "\n")
    return(NULL)
  })
}

cat("\n=== Processing completed for all species ===\n")
