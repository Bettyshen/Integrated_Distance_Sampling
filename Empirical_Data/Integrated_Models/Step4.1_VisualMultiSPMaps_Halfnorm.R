#---
#title: "Visualize Multi-species Prediction Maps (Restricted + non-restircted) with capping values for continuous variables"
#author: "Betty Shen"
#date: "2025-11-18" 
#output: html_document
#---
#Goal: the goal of the script is to predict density of multiple speciesin 2017 using half normal IDS model with capping values for continuous variables and visualize the results
# We will start visualize restricted species first, then visualize species without restricted area
# We will use the grey background for restricted species within Oregon
# This script is aimed for fix distance 500 m as truncation distance

.libPaths(c("/nfs/stak/users/shenf/.conda/envs/r_env_mapping/lib/R/library"))
#Load library
library(unmarked)
library(dplyr)
library(raster)
library(sf)
library(terra)

# ========= Restrive species list (full + restricted species) ========= #
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

# Get restricted species list
raster_restrictedTIFF_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/RestrictedArea/TiffFile"
species_files <- list.files(raster_restrictedTIFF_dir, pattern = "_PredArea.tif$")
restricted_species <- gsub("_PredArea.tif", "", species_files)
cat("Found", length(restricted_species), "species to process:\n")
print(restricted_species)

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
# ============= Load raster data first ============= #
# =============================================================== #

# Load raster
raster <- readRDS("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Raster/xy2017_EnvironmentalRaster_Nov5.rds")
print(names(raster))

# Set target CRS to Oregon State Plane EPSG: 2992
target_crs <- "EPSG:2992"
cat("Target CRS:", target_crs, "\n")

# Check if the raster coordinates are already in the correct CRS
# If not, we need to convert them from WGS84 to Oregon State Plane
cat("Original raster coordinates range:\n")
cat("X range:", range(raster$x), "\n")
cat("Y range:", range(raster$y), "\n")

# If coordinates are in WGS84 (longitude/latitude), convert to Oregon State Plane
# WGS84 coordinates typically range from -180 to 180 for longitude and -90 to 90 for latitude
if (min(raster$x) < -180 || max(raster$x) > 180 || min(raster$y) < -90 || max(raster$y) > 90) {
  cat("Coordinates appear to be in WGS84. Converting to Oregon State Plane...\n")
  
  # Create a temporary raster to convert coordinates
  temp_raster <- rast(raster[, c("x", "y", "EVI_median")], type = "xyz")
  crs(temp_raster) <- "EPSG:4326"  # Set WGS84 CRS
  
  # Convert to Oregon State Plane
  temp_raster_oregon <- project(temp_raster, target_crs)
  
  # Extract the converted coordinates
  raster_coords <- as.data.frame(temp_raster_oregon, xy = TRUE)
  raster$x <- raster_coords$x
  raster$y <- raster_coords$y
  
  cat("Converted coordinates range:\n")
  cat("X range:", range(raster$x), "\n")
  cat("Y range:", range(raster$y), "\n")
  
  # Clean up temporary objects
  rm(temp_raster, temp_raster_oregon, raster_coords)
  gc()
} else {
  cat("Coordinates appear to already be in projected coordinates.\n")
}

print("Load Multi-species IDS model")


# Path to prediction raster (half normal detection models)
raster_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/PredictionMaps/ContinuVars/FixDistance/Halfnorm"

cols <- rev(hcl.colors(100, "Spectral"))



# Read water body raster
water <- rast("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Data/Oregon_Water_Mask.tif")
unique(values(water)) # water is 1, non-water is 0

# =============================================================== #
# ============= For species with restricted area: ============= #
# =============================================================== #

# Process one restricted species at a time to avoid memory issues
# We will use the capping values for each species to visualize the results
# Save visulization raster for each restricted species at the end
for (sp in restricted_species) {
  cat("Processing", sp, "...\n")
  
  # Skip species with NA cap values (insufficient data)
  if (sp %in% species_to_skip) {
    cat("Skipping", sp, "- insufficient data (species has NA cap value)\n")
    next
  }
  
  # Check if species has a cap value
  has_cap_value <- sp %in% names(species_caps)
  if (has_cap_value) {
    cat("Species", sp, "has cap value - will apply capping\n")
  } else {
    cat("Species", sp, "has no cap value - will use raw predictions without capping\n")
  }
  
  sp_filename <- species_to_filename(sp)
  # Load prediction data for this species only - use FullRaster files
  #pred_file <- file.path(raster_dir, paste0("FullRaster_1km_Preds_", sp_filename, ".rds"))
  pred_file <- file.path(raster_dir, paste0("RestrictedRaster_1km_Preds_", sp_filename, ".rds"))

  preds <- readRDS(pred_file)
  
  # Apply capping only if species has a cap value
  preds_capped <- preds
  if (has_cap_value) {
    cap_value <- species_caps[sp]
    print(paste("Capping", sp, "at", cap_value))
    # Cap values for visualization using species-specific cap
    # Only cap the density/prediction values, not coordinates
    # Check what columns are in preds and cap only the prediction column
    if("density_pred" %in% names(preds_capped)) {
      preds_capped$density_pred[preds_capped$density_pred > cap_value] <- cap_value
    } else if("Predicted" %in% names(preds_capped)) {
      preds_capped$Predicted[preds_capped$Predicted > cap_value] <- cap_value
    } else {
      # If no standard prediction column name, cap the first numeric column that's not x or y
      numeric_cols <- names(preds_capped)[sapply(preds_capped, is.numeric)]
      prediction_col <- numeric_cols[!numeric_cols %in% c("x", "y")][1]
      if(!is.na(prediction_col)) {
        preds_capped[[prediction_col]][preds_capped[[prediction_col]] > cap_value] <- cap_value
      }
    }
  } else {
    print(paste("Using raw predictions for", sp, "(no capping applied)"))
  }
  
  # Check the structure of preds data
  cat("Preds data structure:\n")
  print(str(preds_capped))
  cat("Preds dimensions:", dim(preds_capped), "\n")
  cat("Preds column names:", names(preds_capped), "\n")
  
  # Since the prediction files already contain x, y coordinates and density_pred,
  # we can use them directly without any row mismatch checks
  cat("Using prediction data directly with coordinates and predictions.\n")
  raster_out <- preds_capped
  head(raster_out)
  print("Converting to raster")
  
  # Add validation to check for single coordinate values before creating raster
  cat("Checking coordinate data before raster creation:\n")
  cat("Number of unique x coordinates:", length(unique(raster_out$x)), "\n")
  cat("Number of unique y coordinates:", length(unique(raster_out$y)), "\n")
  cat("X range:", range(raster_out$x), "\n")
  cat("Y range:", range(raster_out$y), "\n")
  
  # Check if we have sufficient spatial variation to create a raster
  if (length(unique(raster_out$x)) < 2 || length(unique(raster_out$y)) < 2) {
    cat("ERROR: Insufficient spatial variation detected!\n")
    cat("Cannot create raster with only", length(unique(raster_out$x)), "unique x coordinates and", length(unique(raster_out$y)), "unique y coordinates.\n")
    cat("This usually indicates the restricted area data is too small or has insufficient spatial coverage.\n")
    cat("Skipping", sp, "due to insufficient spatial data.\n")
    next
  }
  
  # Convert to a new raster
  pred_rast <- rast(raster_out[, c("x", "y", "density_pred")], type = "xyz")

  # Check if coordinates are in WGS84 or already projected
  cat("Prediction coordinates range:\n")
  cat("X range:", range(raster_out$x), "\n")
  cat("Y range:", range(raster_out$y), "\n")
  
  # If coordinates look like WGS84 (longitude/latitude), set CRS first then project
  if (min(raster_out$x) >= -180 && max(raster_out$x) <= 180 && 
      min(raster_out$y) >= -90 && max(raster_out$y) <= 90) {
    cat("Coordinates appear to be in WGS84. Setting CRS and projecting...\n")
    crs(pred_rast) <- "EPSG:4326"
    
    # Add error handling for projection
    tryCatch({
      pred_rast <- project(pred_rast, target_crs)
    }, error = function(e) {
      cat("Error during prediction raster projection:", e$message, "\n")
      cat("Setting CRS without projection...\n")
      crs(pred_rast) <<- target_crs
    })
  } else {
    cat("Coordinates appear to be already projected. Setting CRS to target CRS.\n")
    crs(pred_rast) <- target_crs
  }
  
  cat("Predicted raster CRS set to:", crs(pred_rast), "\n")


# ===============Plot and save PNG======================= #
  print("Plotting and saving PNG")
  # List the file path to save the PNG
  png_file <- file.path("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Raster/PredictionMap/ContinVars/FixDistance/Halfnorm", 
                        paste0("RestrictedRaster_1km_Preds_", sp_filename, ".png"))
  
  # Start PNG device
  png(png_file, width=1300, height=900)
  # Increase right margin for legend (default c(4, 4, 2, 2))
    # You can tweak the last value (e.g., 6–10) depending on how wide your legend is.
    op <- par(mar = c(4, 4, 2, 12))  

  # 1. Load the species-specific restricted area file.
  restricted_area_file <- file.path(raster_restrictedTIFF_dir, paste0(sp, "_PredArea.tif"))
  restricted_area_rast <- rast(restricted_area_file)

  # --- FIX for "CRS NOT SET" ERROR ---
  # If the loaded .tif file has no CRS, we must set it before projecting.
  # We'll assume WGS84 (EPSG:4326), a common standard for geographic data.
  if (crs(restricted_area_rast) == "" || is.na(crs(restricted_area_rast))) {
    cat("Warning: CRS for", basename(restricted_area_file), "was not set. Assuming WGS84.\n")
    crs(restricted_area_rast) <- "EPSG:4326"
  }

  # 2. Ensure the mask's CRS and resolution match the prediction raster.
  # Add error handling for projection
  tryCatch({
    # Check if the CRS is valid before projecting
    if (!is.na(crs(restricted_area_rast)) && crs(restricted_area_rast) != "") {
      restricted_area_rast <- project(restricted_area_rast, crs(pred_rast))
      restricted_area_rast <- resample(restricted_area_rast, pred_rast, method="near")
    } else {
      cat("Warning: Invalid CRS for restricted area raster. Skipping projection.\n")
      # If projection fails, try to resample without projection
      restricted_area_rast <- resample(restricted_area_rast, pred_rast, method="near")
    }
  }, error = function(e) {
    cat("Error during projection:", e$message, "\n")
    cat("Attempting to resample without projection...\n")
    # If projection fails, try to resample without projection
    restricted_area_rast <<- resample(restricted_area_rast, pred_rast, method="near")
  })

  # 3. Apply the mask to the prediction data.
  pred_rast_masked <- mask(pred_rast, restricted_area_rast)

  # --- Mask out water bodies ---
 # ============================================================
# Goal: Mask out all waterbodies (with slight buffer) from prediction rasters
# ============================================================


# 2️⃣ Ensure it’s numeric (not logical)
# Use ifel() so output is numeric and keeps CRS/extent properly
water_mask <- ifel(water == 1, 1, 0)

# 3️⃣ Check CRS — fix if missing
if (crs(water_mask) == "" || is.na(crs(water_mask))) {
  cat("Warning: CRS for water mask was not set. Assuming WGS84.\n")
  crs(water_mask) <- "EPSG:4326"
}

# 4️⃣ Match projection and resolution to your prediction raster
#    (replace with the actual prediction raster you're masking)
target_crs <- crs(pred_rast_masked)

water_rast <- project(water_mask, target_crs)
water_rast <- resample(water_rast, pred_rast_masked, method = "near")

# 5️⃣ Buffer (dilate) the water mask to eliminate edge artifacts
#    This expands water pixels outward by N pixels (1 = 1 cell width)
#    You can increase `dilate = 2` for stronger buffering
#buffered_water <- focal(water_rast, w = matrix(1, 3, 3), fun = max, na.policy = "omit")

# Temporarily upscale resolution (half cell size)
fine <- resample(water_rast, disagg(water_rast, 2), method="near")
fine_buf <- focal(fine, w = matrix(1, 3, 3), fun = max)
buffered_water <- aggregate(fine_buf, 2, fun = "max")


# 6️⃣ Apply the mask — remove (set to NA) all buffered waterbody pixels
pred_rast_masked <- mask(pred_rast_masked, buffered_water, maskvalues = 1)


  
  # --- LAYERED PLOTTING ---

  # A. Create the full Oregon grey background
  oregon_raster <- rast("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Raster/EnvironmentalRaster_Nov5.tif")
  
  # Add the same defensive CRS check for the Oregon background raster
  if (crs(oregon_raster) == "" || is.na(crs(oregon_raster))) {
    cat("Warning: CRS for Oregon background raster was not set. Assuming WGS84.\n")
    crs(oregon_raster) <- "EPSG:4326"
  }
  
  # Project the background raster if its CRS doesn't match the target CRS
  # Add error handling for projection
  tryCatch({
    if(crs(oregon_raster) != target_crs && !is.na(crs(oregon_raster)) && crs(oregon_raster) != "") {
      oregon_raster <- project(oregon_raster, target_crs)
    }
  }, error = function(e) {
    cat("Error during Oregon background projection:", e$message, "\n")
    cat("Continuing with original CRS...\n")
  })
  oregon_mask <- oregon_raster[[1]]
  oregon_mask[!is.na(oregon_mask)] <- 1

  # B. Plot the grey Oregon background layer first
  plot(oregon_mask,
       col = 'lightgrey',
       legend = FALSE,
       main = filename_to_species(sp_filename),
       axes = TRUE)

  # C. Add the correctly masked species data on top
  plot(pred_rast_masked,
       col = cols,
       add = TRUE,
       legend = TRUE)
  # Restore margins
  par(op) # Restore original margins
  # Close the PNG device
  dev.off()
  
  # Clear memory after processing each species
  rm(preds, raster_out, pred_rast, restricted_area_rast, pred_rast_masked, oregon_raster, oregon_mask)
  gc()
  
  cat("Done processing", sp, "\n")
}


# =============================================================== #
# ============= For other species (non-restricted): ============= #
# =============================================================== #

# List of selected species (all species except restricted species)
species <- All_species[!All_species %in% restricted_species]
cat("Found", length(species), "species to process:\n")
print(species)


# Process and visualize each species in one efficient loop
# Process and visualize each species in one efficient loop
for (sp in species) {
  cat("Processing and visualizing", sp, "...\n")
  
  # Check if prediction file exists
  pred_file <- file.path(raster_dir, paste0("FullRaster_1km_Preds_", gsub(' ', '_', sp), ".rds"))
  
  if (!file.exists(pred_file)) {
    cat("Warning: Prediction file not found for", sp, "\n")
    cat("File path:", pred_file, "\n")
    cat("Skipping", sp, "and moving to next species...\n")

    next
  }
  
  # Load prediction data for this species only
  preds <- readRDS(pred_file)
  
  # Check the structure of preds to determine how to handle it
  cat("Preds data structure:\n")
  print(str(preds))
  cat("Preds column names:", names(preds), "\n")
  
  # If preds already has x, y, and density_pred columns, use it directly
  # Otherwise, merge with raster (for older format files)
  if ("x" %in% names(preds) && "y" %in% names(preds) && "density_pred" %in% names(preds)) {
    raster_out <- preds
    pred_col <- "density_pred"
  } else if ("x" %in% names(preds) && "y" %in% names(preds) && "Predicted" %in% names(preds)) {
    raster_out <- preds
    pred_col <- "Predicted"
  } else {
    # Old format: merge with raster
    raster_out <- cbind(raster[, c("x", "y")], preds)
    # Determine which prediction column exists
    if ("density_pred" %in% names(raster_out)) {
      pred_col <- "density_pred"
    } else if ("Predicted" %in% names(raster_out)) {
      pred_col <- "Predicted"
    } else {
      # Find first numeric column that's not x or y
      numeric_cols <- names(raster_out)[sapply(raster_out, is.numeric)]
      pred_col <- numeric_cols[!numeric_cols %in% c("x", "y")][1]
      if (is.na(pred_col)) {
        stop("Could not find prediction column in data")
      }
    }
  }
  
  head(raster_out)
  print("Converting to raster")
  # Convert to a new raster
  pred_rast <- rast(raster_out[, c("x", "y", pred_col)], type = "xyz")
  
  # Check if coordinates are in WGS84 or already projected
  cat("Prediction coordinates range:\n")
  cat("X range:", range(raster_out$x), "\n")
  cat("Y range:", range(raster_out$y), "\n")
  
  # If coordinates look like WGS84 (longitude/latitude), set CRS first then project
  if (min(raster_out$x) >= -180 && max(raster_out$x) <= 180 && 
      min(raster_out$y) >= -90 && max(raster_out$y) <= 90) {
    cat("Coordinates appear to be in WGS84. Setting CRS and projecting...\n")
    crs(pred_rast) <- "EPSG:4326"
    
    # Add error handling for projection
    tryCatch({
      pred_rast <- project(pred_rast, target_crs)
    }, error = function(e) {
      cat("Error during prediction raster projection:", e$message, "\n")
      cat("Setting CRS without projection...\n")
      crs(pred_rast) <<- target_crs
    })
  } else {
    cat("Coordinates appear to be already projected. Setting CRS to target CRS.\n")
    crs(pred_rast) <- target_crs
  }
  
  cat("Predicted raster CRS set to:", crs(pred_rast), "\n")

  # ===============Plot and save PNG======================= #
  print("Plotting and saving PNG")
  # List the file path to save the PNG
  png_file <- file.path("/nfs/stak/users/shenf/hpc-share/IntegrationModel/Raster/PredictionMap/ContinVars/FixDistance/Halfnorm", 
                        paste0("FullRaster_1km_Preds_", gsub(' ', '_', sp), ".png"))
  
  # For non-restricted species, plot the full prediction raster without Oregon boundary mask
  # since they should be predicted throughout the entire state
  png(png_file, width=1200, height=900)
  plot(pred_rast, main = paste(sp), col = cols)
  dev.off()
  
  # Clear memory after processing each species
  #rm(preds, raster_out, pred_rast, oregon_raster)
  rm(preds, raster_out, pred_rast)
  gc()
  
  cat("Done processing and visualizing", sp, "\n")
}
