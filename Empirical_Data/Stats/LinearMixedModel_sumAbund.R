# Goal: Test if half-normal and species-specific detection functions predicted sum abundance are significantly different
# Date: 2026-01-16
# Author: Betty Shen
# 
# KEY MODEL STRUCTURE CHANGE:
# Uses a SINGLE shared observation error term (Obs.Err) instead of species-specific errors.
# This is critical because:
# 1. We have only 2 observations per species (one for each model type)
# 2. These are deterministic model outputs, not noisy observations
# 3. With species-specific error and only 2 observations, the error can absorb the signal,
#    making Beta1 (the difference between models) appear non-significant
# Remove environmnet
rm(list = ls())
# Load required libraries
library(nimble)
library(tidyverse)
library(dplyr)
library(coda)
library(data.table)
library(mcmcplots)
library(MCMCvis)
source('/Users/shenf/Documents/Bayes_Class/Week2_Nimble/attach.nimble_v2.R')

# Load data
data <- read.csv("/Volumes/T7/Chapter2/DataIntegration/Results/Abundance_Compare/Abundance_Statistics_Comparison_merged.csv")
str(data)

# Remove Mountain Quail data and best key function is not half-normal
data.select <- data[data$Common_Name != "Mountain Quail" & data$Best_Key_Function != "half_normal", ]
nrow(data.select)
head(data.select, 10)
names(data.select)
# Select columns of interest
data.cols <- data.select %>%
    select(Common_Name, HalfNorm_Sum, SpeciesSpec_Sum)
str(data.cols)

# Convert data into long format: species, model type (half-normal or species-specific), and abundance
data_long <- data.cols %>%
    pivot_longer(
        cols = c(HalfNorm_Sum, SpeciesSpec_Sum),
        names_to = "Model_Type",
        values_to = "Abundance"
    ) %>%
    mutate(
        Model_Type = case_when(
            Model_Type == "HalfNorm_Sum" ~ "Half-normal",
            Model_Type == "SpeciesSpec_Sum" ~ "Species-specific",
            TRUE ~ Model_Type
        )
    )
str(data_long)
head(data_long, 10)

# Compose linear mixed-effect model - by creating a mean effect that the species-level effects are distributed around it

model <- nimbleCode({
    # Priors
    Int.hat ~ dt(mu = 0, sigma = 1, df = 1)
    Beta1.hat ~ dt(mu = 0, sigma = 1, df = 1)
    Sigma.int ~ T(dt(mu = 0, sigma = 1, df = 1), 0,) # Half Cauchy
    Sigma.beta1 ~ T(dt(mu = 0, sigma = 1, df = 1), 0,) # Half Cauchy
    # Use a SINGLE observation error term (not species-specific)
    # This is appropriate because we're comparing deterministic model outputs
    # and with only 2 observations per species, species-specific error can absorb the signal
    Obs.Err ~ T(dt(mu = 0, sigma = 1, df = 1), 0,) # Half Cauchy

    for (j in 1:n.species) {
        Intercept[j] ~ dnorm(mean = Int.hat, sd = Sigma.int)
        Beta1[j] ~ dnorm(mean = Beta1.hat, sd = Sigma.beta1)
    } # end j loop

for (i in 1:n.obs) {


    # Process model
    Abundance.exp[i] <- Intercept[Species[i]] + Beta1[Species[i]]*Model.Type[i]

    # Observation model (likelihood)
    # Use shared observation error across all species
    Abundance[i] ~ dnorm(mean = Abundance.exp[i], sd = Obs.Err)

} # end i loop
   
    
})

# Parameters to monitor (including hyperparameters for better diagnostics)
# Note: Obs.Err is now a scalar, not a vector
parameters <- c("Intercept", "Beta1", "Obs.Err", "Int.hat", "Beta1.hat", "Sigma.int", "Sigma.beta1")

# MCMC settings - INCREASED for better convergence
ni <- 500000  # Increased from 100000
nt <- 40      
nb <- 250000  # Increased burn-in from 50000
nc <- 3

# Convert Model_Type to numeric: Half-normal = 0, Species-specific = 1
data_long$Model_Type_numeric <- ifelse(data_long$Model_Type == "Half-normal", 0, 1)
print(data_long)

# Check for missing values and scale issues
cat("\n=== DATA DIAGNOSTICS ===\n")
cat("Number of observations:", nrow(data_long), "\n")
cat("Number of species:", length(unique(data_long$Common_Name)), "\n")
cat("Abundance range:", range(data_long$Abundance, na.rm = TRUE), "\n")
cat("Abundance mean:", mean(data_long$Abundance, na.rm = TRUE), "\n")
cat("Abundance SD:", sd(data_long$Abundance, na.rm = TRUE), "\n")
cat("Missing values:", sum(is.na(data_long$Abundance)), "\n")

# Remove any missing values if present
if (sum(is.na(data_long$Abundance)) > 0) {
    cat("Warning: Removing", sum(is.na(data_long$Abundance)), "rows with missing abundance values\n")
    data_long <- data_long[!is.na(data_long$Abundance), ]
}

# Data 
nimble.data <- list(Model.Type = data_long$Model_Type_numeric,
                    Abundance = c(scale(data_long$Abundance)))
print(nimble.data)
# Create factor for species (store this to ensure consistent ordering)
species_factor <- factor(data_long$Common_Name)
species_names_ordered <- levels(species_factor)  # Store ordered species names

nimble.constants <- list(n.obs = nrow(data_long),
                         n.species = length(unique(data_long$Common_Name)),
                         Species = as.numeric(species_factor))

# Create initial values function for better starting points
# This helps chains start from different, reasonable values
# NOTE: Since abundance data is scaled (mean ≈ 0, SD ≈ 1), we use scaled statistics
inits <- function() {
    # Get mean from SCALED data (should be ≈ 0)
    abundance_mean <- mean(nimble.data$Abundance, na.rm = TRUE)
    
    # Initialize hyperparameters using scaled data statistics
    # Since data is scaled, intercept should be near 0, and SDs should be on scale of 1
    Int.hat.init <- rnorm(1, mean = abundance_mean, sd = 0.5)  # Near 0 for scaled data
    Beta1.hat.init <- rnorm(1, mean = 0, sd = 0.2)              # Small effect for scaled data
    Sigma.int.init <- abs(rnorm(1, mean = 0.5, sd = 0.2))      # Variance on scale of 1
    Sigma.beta1.init <- abs(rnorm(1, mean = 0.2, sd = 0.1))    # Smaller variance for slope
    
    # Initialize species-level parameters
    Intercept.init <- rnorm(nimble.constants$n.species, 
                           mean = Int.hat.init, 
                           sd = Sigma.int.init * 0.5)
    Beta1.init <- rnorm(nimble.constants$n.species,
                       mean = Beta1.hat.init,
                       sd = Sigma.beta1.init * 0.5)
    # Observation error should be on scale of scaled data (SD ≈ 1)
    # Now a single scalar value (not species-specific)
    Obs.Err.init <- abs(rnorm(1,
                             mean = 0.5,  # Reasonable for scaled data
                             sd = 0.2))
    
    list(Int.hat = Int.hat.init,
         Beta1.hat = Beta1.hat.init,
         Sigma.int = Sigma.int.init,
         Sigma.beta1 = Sigma.beta1.init,
         Intercept = Intercept.init,
         Beta1 = Beta1.init,
         Obs.Err = Obs.Err.init)
}

cat("\n=== RUNNING MCMC (this may take a while) ===\n")
cat("Iterations:", ni, "\n")
cat("Burn-in:", nb, "\n")
cat("Thinning:", nt, "\n")
cat("Chains:", nc, "\n")
cat("Effective samples per chain:", (ni - nb) / nt, "\n\n")

# Run MCMC with adaptive settings
mcmc.output <- nimbleMCMC(code = model,
                        data = nimble.data,
                        constants = nimble.constants,
                        inits = inits,
                        monitors = parameters,
                        nchains = nc,
                        niter = ni,
                        nburnin = nb,
                        thin = nt,
                        summary = TRUE,
                        samplesAsCodaMCMC = TRUE,
                        WAIC = TRUE,
                        # Adaptive MCMC settings for better mixing
                        setSeed = TRUE,
                        progressBar = TRUE)

attach.nimble(mcmc.output$samples)

MCMCtrace(object = mcmc.output$samples,
    pdf = FALSE,
    ind = TRUE,
    params = c("Intercept", "Beta1", "Obs.Err"))
# Note: Obs.Err is now a scalar parameter, not a vector

# ===== CONVERGENCE DIAGNOSTICS =====
cat("\n=== CONVERGENCE DIAGNOSTICS ===\n")

# Gelman-Rubin diagnostic (Rhat or PSRF)
# Values should be < 1.1 for convergence, < 1.05 is ideal
gelman.diag.result <- gelman.diag(mcmc.output$samples)
print(gelman.diag.result) # 1.04

# Check which parameters have Rhat > 1.1
rhat.values <- gelman.diag.result$psrf[, "Point est."]
cat("\nParameters with Rhat > 1.1 (poor convergence):\n")
print(names(rhat.values[rhat.values > 1.1]))
cat("\nParameters with Rhat > 1.05 (marginal convergence):\n")
print(names(rhat.values[rhat.values > 1.05 & rhat.values <= 1.1]))

# Effective sample size
cat("\n=== EFFECTIVE SAMPLE SIZE ===\n")
effectiveSize.result <- effectiveSize(mcmc.output$samples)
print(effectiveSize.result)
cat("\nMinimum effective sample size:", min(effectiveSize.result), "\n")
cat("Parameters with ESS < 400 (may need more samples):\n")
print(names(effectiveSize.result[effectiveSize.result < 400]))

# Geweke diagnostic (tests for convergence within chains)
cat("\n=== GEWEKE DIAGNOSTIC (within-chain convergence) ===\n")
geweke.result <- geweke.diag(mcmc.output$samples)
print(geweke.result)

# Visualize all of the relevant plots at the same time
mcmcplot(mcmc.output$samples)

# Save MCMC output
saveRDS(mcmc.output, file = "/Volumes/T7/Chapter2/DataIntegration/Results/Bayes/LinearMixedModel_SumAbund_output.rds")

# Read MCMC output
mcmc.output <- readRDS("/Volumes/T7/Chapter2/DataIntegration/Results/Bayes/LinearMixedModel_SumAbund_output.rds")
attach.nimble(mcmc.output$samples)
# ===== EXTRACT SPECIES-SPECIFIC BETA1 (EXCLUDE HYPERPARAMETERS) =====
# Beta1 matrix includes species-specific values AND hyperparameters (Beta1.hat, Sigma.beta1)
# We only want the species-specific values (first n.species columns)

cat("\n=== EXTRACTING SPECIES-SPECIFIC BETA1 ===\n")
cat("Total Beta1 columns (including hyperparameters):", ncol(Beta1), "\n")

# Get number of species from constants (or calculate from data if constants not available)
if (exists("nimble.constants")) {
    n_species <- nimble.constants$n.species
} else {
    n_species <- length(unique(data_long$Common_Name))
}

cat("Number of species:", n_species, "\n")
cat("Extracting first", n_species, "columns (species-specific Beta1)\n")
cat("Ignoring hyperparameters (Beta1.hat, Sigma.beta1)\n\n")

# Extract only species-specific Beta1 columns (first n_species columns)
Beta1_species <- Beta1[, 1:n_species, drop = FALSE]

# Summarize Posterior Distribution for species-specific Beta1
beta1_summary <- apply(Beta1_species, 2, quantile, probs = c(0.025, 0.5, 0.975))
print(beta1_summary)

# 95% credible interval for Beta1
beta1_ci <- apply(Beta1_species, 2, quantile, probs = c(0.025, 0.975))
print(beta1_ci)

# ===== CHECK IF BETA1 95% CRI OVERLAPS 0 =====
# Get species names in the same order as Beta1 columns
# Use the stored species_names_ordered if available (from when constants were created)
# Otherwise, recreate from the factor

cat("\n=== SPECIES NAME EXTRACTION ===\n")

if (exists("species_names_ordered") && length(species_names_ordered) == n_species) {
    # Use the stored species names (best option - matches model exactly)
    species_names <- species_names_ordered
    cat("Using stored species names from model constants:", length(species_names), "species\n")
} else {
    # Fallback: recreate from data
    cat("Recreating species names from data...\n")
    species_factor_temp <- factor(data_long$Common_Name)
    species_names <- levels(species_factor_temp)
    cat("Number of unique species from data:", length(species_names), "\n")
    
    if (length(species_names) != n_species) {
        # Try to get species names in the exact order they appear in the factor
        species_names_ordered <- unique(data_long$Common_Name)[order(unique(as.numeric(species_factor_temp)))]
        
        if (length(species_names_ordered) == n_species) {
            species_names <- species_names_ordered
            cat("Fixed: Using ordered unique species names\n")
        } else {
            # If still mismatched, this suggests data/model mismatch
            stop(paste("Cannot match species names to Beta1 columns.\n",
                      "Species count:", length(species_names), 
                      "Expected species:", n_species,
                      "\nThis may indicate the model was run with different data."))
        }
    }
}

cat("Final number of species names:", length(species_names), "\n\n")

# Final safety check: ensure lengths match exactly
if (length(species_names) != ncol(Beta1_species)) {
    stop(paste("Cannot proceed: species_names length (", length(species_names), 
               ") does not match Beta1_species columns (", ncol(Beta1_species), ")"))
}

# Create a data frame with Beta1 results
beta1_results <- data.frame(
    Species = species_names,
    Lower_2.5 = beta1_ci[1, ],
    Upper_97.5 = beta1_ci[2, ],
    Sum = beta1_summary[2, ],
    Overlaps_Zero = (beta1_ci[1, ] < 0 & beta1_ci[2, ] > 0),
    Significant = !(beta1_ci[1, ] < 0 & beta1_ci[2, ] > 0)  # Significant if CRI does NOT overlap 0
)

# Sort by median Beta1 for easier interpretation
beta1_results <- beta1_results[order(beta1_results$Sum), ]

cat("\n=== BETA1 RESULTS: 95% CREDIBLE INTERVALS ===\n")
cat("Note: 'Overlaps_Zero' = TRUE means the 95% CRI includes 0 (effect not significant)\n")
cat("      'Significant' = TRUE means the 95% CRI excludes 0 (effect is significant)\n\n")
print(beta1_results)

# Summary statistics
cat("\n=== SUMMARY ===\n")
cat("Number of species with Beta1 95% CRI overlapping 0:", sum(beta1_results$Overlaps_Zero), "\n")
cat("Number of species with Beta1 95% CRI NOT overlapping 0 (significant):", sum(beta1_results$Significant), "\n")
cat("Species with significant positive effect (lower bound > 0):", 
    sum(beta1_results$Lower_2.5 > 0), "\n")
cat("Species with significant negative effect (upper bound < 0):", 
    sum(beta1_results$Upper_97.5 < 0), "\n")

# Visualize which species have significant effects
cat("\n=== SPECIES WITH SIGNIFICANT EFFECTS (95% CRI excludes 0) ===\n")
significant_species <- beta1_results[beta1_results$Significant, ]
if (nrow(significant_species) > 0) {
    print(significant_species[, c("Species", "Lower_2.5", "Sum", "Upper_97.5")])
} else {
    cat("No species have significant effects (all 95% CRIs overlap 0)\n")
}
# ===== CREATE EFFECT PLOTS FOR BETA1 =====
library(ggplot2)

# Read taxonomy checklist
clements <- read.csv("/Volumes/T7/Chapter2/DataIntegration/Data/Taxonomy/Clements_v2025.csv")
cat("\n=== CLEMENTS CHECKLIST INFO ===\n")
cat("Columns in Clements data:", paste(names(clements), collapse = ", "), "\n")
cat("Number of species in Clements:", nrow(clements), "\n")

# Check if English.name column exists
if (!"English.name" %in% names(clements)) {
    stop("Error: 'English.name' column not found in Clements data. Available columns: ", 
         paste(names(clements), collapse = ", "))
}

# Prepare data for plotting - Order species according to Clements checklist
cat("\n=== ORDERING SPECIES BY CLEMENTS CHECKLIST ===\n")
cat("Number of species in beta1_results:", nrow(beta1_results), "\n")

# Match species names from beta1_results with Clements checklist
# Create a mapping: get the order index from Clements for each species
species_order <- match(beta1_results$Species, clements$English.name)

# Check for species not found in Clements
missing_species <- beta1_results$Species[is.na(species_order)]
if (length(missing_species) > 0) {
    cat("Warning: The following species were not found in Clements checklist:\n")
    print(missing_species)
    cat("\nThese will be placed at the end in their original order.\n")
}

# For species found in Clements, order by their position in the checklist
# For species not found, assign a large number so they appear at the end
species_order[is.na(species_order)] <- max(species_order, na.rm = TRUE) + 1:sum(is.na(species_order))

# Sort beta1_results by Clements order
beta1_results <- beta1_results[order(species_order), ]

cat("Species ordered by Clements checklist.\n")
cat("First 5 species:", paste(beta1_results$Species[1:5], collapse = ", "), "\n")
cat("Last 5 species:", paste(beta1_results$Species[(nrow(beta1_results)-4):nrow(beta1_results)], collapse = ", "), "\n")

# Create factor with species in Clements order (reversed for plotting: top to bottom)
# This ensures species appear from top to bottom in Clements order
beta1_results$Species <- factor(beta1_results$Species, levels = rev(beta1_results$Species))

# Split species into two groups for two figures
n_species_total <- nrow(beta1_results)
split_point <- ceiling(n_species_total / 2)

# Group 1: First half of species (maintain Clements order)
beta1_plot1 <- beta1_results[1:split_point, ]
# Preserve factor levels to maintain Clements order in plot
beta1_plot1$Species <- factor(beta1_plot1$Species, 
                             levels = levels(beta1_results$Species)[levels(beta1_results$Species) %in% beta1_plot1$Species])

# Group 2: Second half of species (maintain Clements order)
beta1_plot2 <- beta1_results[(split_point + 1):n_species_total, ]
# Preserve factor levels to maintain Clements order in plot
beta1_plot2$Species <- factor(beta1_plot2$Species, 
                              levels = levels(beta1_results$Species)[levels(beta1_results$Species) %in% beta1_plot2$Species])

# Function to create effect plot
create_effect_plot <- function(data, title_suffix) {
    p <- ggplot(data, aes(x = Sum, y = Species)) +
        # Add vertical line at 0 for reference
        geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
        # Add 95% CRI as error bars
        geom_errorbarh(aes(xmin = Lower_2.5, xmax = Upper_97.5), 
                      height = 0.3, linewidth = 0.5, color = "gray40") +
        # Add median point
        geom_point(aes(x = Sum), size = 2, color = "orange", fill = "orange", 
                  shape = 21, stroke = 1) +

        labs(
            x = "Estimated Coefficient",
            y = "Species",
            #title = paste("Beta1 Effect Estimates with 95% Credible Intervals", title_suffix)
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold"),
            panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            legend.position = "bottom"
        )
    return(p)
}

# Create the two plots
cat("\n=== CREATING EFFECT PLOTS ===\n")
plot1 <- create_effect_plot(beta1_plot1, "(Species 1-33)")
plot2 <- create_effect_plot(beta1_plot2, paste0("(Species ", split_point + 1, "-", n_species_total, ")"))

# Display plots
print(plot1)
print(plot2)

# Save plots
plot_dir <- "/Volumes/T7/Chapter2/DataIntegration/Results/Bayes/"
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
}

ggsave(filename = paste0(plot_dir, "Beta1_SumAbundEffectPlot_Part1.png"), 
       plot = plot1, width = 10, height = 12, dpi = 300)
cat("Saved: Beta1_SumAbundEffectPlot_Part1.png\n")

ggsave(filename = paste0(plot_dir, "Beta1_SumAbundEffectPlot_Part2.png"), 
       plot = plot2, width = 10, height = 12, dpi = 300)
cat("Saved: Beta1_SumAbundEffectPlot_Part2.png\n")
cat("Plots saved to:", plot_dir, "\n")
