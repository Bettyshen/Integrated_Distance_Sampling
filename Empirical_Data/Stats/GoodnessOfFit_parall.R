#*This script will be used on cluster*
  
#Goal: Calculate goodness-of-fitdifferences by using half-normal vs. species-specific detection key function via eBird + Oregon 2020 from zero-filled data for multiple species

#- Extract goodness-of-fit measures for predictors

#=== Load library & Data ====#
library(unmarked)
library(dplyr)
library(tidyr)
library(AICcmodavg)

# Base directory
base_dir <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Data"

# Oregon 2020 
OR2020 <- read.csv(sprintf("%s/Oregon2020_BirdEnvs_zf_Trunc.csv", base_dir))

# eBird
eBird <- read.csv(sprintf("%s/eBird_zf_envr.csv", base_dir))

# Key function
key <- read.csv(sprintf("%s/Best_key_functions_trunc.csv", base_dir))

#======== Scale environmental covariates ==========#
# List predictors we want to scale
predictors <- c("EVI_median","SR_B1_median","SR_B5_median","ppt_median", "tmean_median")

# Add scaled versions of those columns with a "_scaled" suffix in new columns
# Oregon 2020
oregon2020_scaled <- OR2020 %>%
  mutate(across(all_of(predictors), 
                ~ as.numeric(scale(.)), 
                .names = "{.col}_scaled"))
names(oregon2020_scaled)
# eBird
eBird_scaled <- eBird %>%
  mutate(across(all_of(predictors),
                ~ as.numeric(scale(.)),
                .names = "{.col}_scaled"))

#======== Balance dataset in advance =================#
# - Detections vs. non-detections
balance_bird_occurrence <- function(dataset, detect_col) {
  set.seed(123)
  balanced_data <- data.frame()  # Initialize empty result
  
  for (sp in unique(dataset$Common_Name)) {
    
    sp.data <- dataset[dataset$Common_Name == sp,]
    
    bird_detect <- sp.data[sp.data[[detect_col]] > 0, ]
    bird_nonDetect <- sp.data[sp.data[[detect_col]] == 0,]
    
    n.detect <- nrow(bird_detect)
    n.nonDetect <- nrow(bird_nonDetect)
    
    # Only balance if there are both detections and non-detections
    if (n.detect > 0 && n.nonDetect > 0) {
      n.sample <- n.detect
      
      nonDetect_sample <- bird_nonDetect[sample(nrow(bird_nonDetect), n.sample), ]
      
      sp_balanced <- rbind(bird_detect, nonDetect_sample)
      
      balanced_data <- rbind(balanced_data, sp_balanced)
    }
  }
  
  return (balanced_data)
  
}

# Balance Oregon 2020 & eBird
OR2020_balanced <- balance_bird_occurrence(oregon2020_scaled, "Occur")
eBird_balanced <- balance_bird_occurrence(eBird_scaled, "observation_count")

#Loop through each species and export goodness-of-fit measures using parallel processing

# All species (unique)
unique_SP <- unique(OR2020_balanced$Common_Name)

# Scaled predictor variables
predictors_scaled <- c("EVI_median_scaled","SR_B1_median_scaled",
                       "SR_B5_median_scaled","ppt_median_scaled", 
                       "tmean_median_scaled")

# Load parallel processing libraries
library(parallel)
library(doParallel)
library(foreach)

# Set up parallel processing
num_cores <- detectCores() - 1  # Leave one core free for system
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Initialize the results CSV file with headers
results_file <- "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Results/Goodness_of_Fit_incremental.csv"
headers <- data.frame(Common_Name = character(),
                     Number_of_bootstrap = integer(),
                     P_value_Pcount = numeric(),
                     P_value_DS_half = numeric(),
                     P_value_DS_Spec = numeric(),
                     Estimate_C_hat_Pcount = numeric(),
                     Estimate_C_hat_DS_half = numeric(),
                     Estimate_C_hat_DS_Spec = numeric())
write.csv(headers, results_file, row.names = FALSE)

# Process species in parallel with incremental saving
all_birds_results <- foreach(sp = unique_SP, .packages = c("unmarked", "dplyr", "tidyr", "AICcmodavg"), 
                            .combine = 'c', .errorhandling = 'pass') %dopar% {
  
  set.seed(123)
  
  # Process species
  print(paste("Processing species:", sp))
  
  tryCatch({
    ##====== Oregon 2020 data =====##
    sp_data <- OR2020_balanced %>% filter(Common_Name == sp)
    unique_sp_data <- sp_data %>%
      group_by(Unique_Observation_ID) %>%
      slice_sample(n = 1) %>%
      ungroup() %>% as.data.frame()
    
    trunc_dist <- key %>% filter(Common_Name == sp) %>% pull(Truncation_Distance)
    if (length(trunc_dist) == 0) {
      Goodness_out <- data.frame(Common_Name = sp,
                                Number_of_bootstrap = NA,
                                P_value_Pcount = NA,
                                P_value_DS_half = NA,
                                P_value_DS_Spec = NA,
                                Estimate_C_hat_Pcount = NA,
                                Estimate_C_hat_DS_half = NA,
                                Estimate_C_hat_DS_Spec = NA
      )
      return(list(Goodness_out = Goodness_out))
    }
    
    unique_sp_data$Unique_Observation_ID <- as.factor(unique_sp_data$Unique_Observation_ID)
    dist_breaks <- seq(0, trunc_dist + 10, by = 15)
    dis_sp <- formatDistData(unique_sp_data, distCol = "Distance",
                             transectNameCol = "Unique_Observation_ID",
                             dist.breaks = dist_breaks)
    
    covs_sp <- unique_sp_data[, c("Unique_Observation_ID", predictors_scaled)]
    umf_sp <- unmarkedFrameDS(y = as.matrix(dis_sp), siteCovs = covs_sp, 
                              survey = "point", dist.breaks = dist_breaks, unitsIn = "m")
    
    ##====== eBird data ======##
    sp_eB <- eBird_balanced %>% filter(Common_Name == sp)
    sp_eB$sampling_event_identifier <- as.factor(sp_eB$sampling_event_identifier)
    covs_eB <- sp_eB[, c("sampling_event_identifier", predictors_scaled)]
    Mpc <- nrow(sp_eB)
    db2 <- c(0, trunc_dist)
    umf_eB <- unmarkedFrameDS(y = matrix(sp_eB$observation_count, Mpc, 1), 
                              siteCovs = covs_eB, survey = "point", dist.breaks = db2, unitsIn = "m")
    umf_pc <- unmarkedFramePCount(y = umf_eB@y, siteCovs = umf_eB@siteCovs)
    
    # === Best key function of the species === #
    key_function <- key %>% filter(Common_Name == sp) %>% pull(Best_Key_Function)
    
    # Fit models
    # ====== Point Count Data (eBird) ====== #
    bird.Pcount <- pcount( ~1 ~ EVI_median_scaled + SR_B1_median_scaled + SR_B5_median_scaled +
                             ppt_median_scaled + tmean_median_scaled,
                          data = umf_pc)

    # ====== Distance Sampling Data (Oregon 2020) - Half-normal key function ====== #
    bird.DS.half <- distsamp(~1 ~EVI_median_scaled + SR_B1_median_scaled + SR_B5_median_scaled +
                             ppt_median_scaled + tmean_median_scaled,
                             data = umf_sp, keyfun = "halfnorm")
    # ====== Distance Sampling Data (Oregon 2020) - Species-specific key function ====== #
    bird.DS.Spec <- distsamp(~1 ~EVI_median_scaled + SR_B1_median_scaled + SR_B5_median_scaled +
                             ppt_median_scaled + tmean_median_scaled,
                             data = umf_sp, keyfun = key_function)

    # Calculate goodness-of-fit measures
    print("Calculating goodness-of-fit measures")
    # Set number of simulations
    nsim <- 100  # You can adjust this value as needed
    
    # ====== Point Count Data (eBird) ====== #
    bird.Pcount.Goodness <- Nmix.gof.test(bird.Pcount, nsim = nsim, parallel = FALSE)  # Set to FALSE since we're already parallelizing at species level
    # ====== Distance Sampling Data (Oregon 2020) - Half-normal key function ====== #
    bird.DS.half.Goodness <- Nmix.gof.test(bird.DS.half, nsim = nsim, parallel = FALSE)
    # ====== Distance Sampling Data (Oregon 2020) - Species-specific key function ====== #
    bird.DS.Spec.Goodness <- Nmix.gof.test(bird.DS.Spec, nsim = nsim, parallel = FALSE)
    
    # Add NA if the model was not trained successfully
    if (is.null(bird.Pcount) || is.null(bird.DS.half) || is.null(bird.DS.Spec)) {
      Goodness_out <- data.frame(Common_Name = sp,
      Number_of_bootstrap = NA,
      P_value_Pcount = NA,
      P_value_DS_half = NA,
      P_value_DS_Spec = NA,
      Estimate_C_hat_Pcount = NA,
      Estimate_C_hat_DS_half = NA,
      Estimate_C_hat_DS_Spec = NA
   )
      
      # Append this species result to the CSV file
      write.table(Goodness_out, results_file, append = TRUE, sep = ",", 
                  col.names = FALSE, row.names = FALSE)
      
      return(list(Goodness_out = Goodness_out))
    }

    
    # Extract goodness-of-fit information both species-specific vs. half-normal key function if trained successfully
    Goodness_out <- data.frame(Common_Name = sp,
                                Number_of_bootstrap = nsim,
                                P_value_Pcount = bird.Pcount.Goodness$p.value,
                                P_value_DS_half = bird.DS.half.Goodness$p.value,
                                P_value_DS_Spec = bird.DS.Spec.Goodness$p.value,
                                Estimate_C_hat_Pcount = bird.Pcount.Goodness$c.hat.est,
                                Estimate_C_hat_DS_half = bird.DS.half.Goodness$c.hat.est,
                                Estimate_C_hat_DS_Spec = bird.DS.Spec.Goodness$c.hat.est
                                )
    
    # Append this species result to the CSV file
    write.table(Goodness_out, results_file, append = TRUE, sep = ",", 
                col.names = FALSE, row.names = FALSE)
    
    return(list(Goodness_out = Goodness_out))
    
  }, error = function(e) {
    # Handle errors gracefully
    print(paste("Error processing species:", sp, "-", e$message))
    Goodness_out <- data.frame(Common_Name = sp,
                              Number_of_bootstrap = NA,
                              P_value_Pcount = NA,
                              P_value_DS_half = NA,
                              P_value_DS_Spec = NA,
                              Estimate_C_hat_Pcount = NA,
                              Estimate_C_hat_DS_half = NA,
                              Estimate_C_hat_DS_Spec = NA
    )
    
    # Append this species result to the CSV file
    write.table(Goodness_out, results_file, append = TRUE, sep = ",", 
                col.names = FALSE, row.names = FALSE)
    
    return(list(Goodness_out = Goodness_out))
  })
}

# Stop parallel cluster
stopCluster(cl)

# Loop completed
print("Parallel processing completed")

# Combine results
Goodness_results <- do.call(rbind, lapply(all_birds_results, function(x) x$Goodness_out))

# Save final combined results
write.csv(Goodness_results, "/nfs/stak/users/shenf/hpc-share/IntegrationModel/Results/Goodness_of_Fit_final.csv", row.names = FALSE)

print("All processing completed successfully!")

