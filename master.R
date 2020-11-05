# Master file to run stock recruit analysis

# -------------------------------------------#
# Load packages, data, and functions
# -------------------------------------------#

setwd("C:/github/fishyTMB") # temporary, set working directory

# read in required libraries
library(TMB)
library(ggplot2)
library(dplyr)
library(tidyr)

# source functions and figure functions
source("functions.r")
source("figures.r")

# download data and save to .csv if not already present, load and format data
source("load_data.r")

# make vector of models to use
models <- c("ricker_basic",
               "ricker_multi_CUs",
               "ricker_SMSY_Sgen",
               "Aggregate_LRPs")

# -------------------------------------------#
# Compile TMB models and load
# -------------------------------------------#

# Compile and load TMB files
for(i in 1:length(models)) {
  compile(paste0("TMB/", models[i], ".cpp"))
  dyn.load(dynlib(paste0("TMB/", models[i])))
}

# -------------------------------------------#
# Prepare data for TMB models
# -------------------------------------------#

# make list of model input lists
model_input_list <- lapply(X=models, make_model_input, SRdat = dat) # model_name argument automatically evaluates to elements of models object
names(model_input_list) <- models # name elements of list after model names

# loadTMB("ricker_basic") # function to compile and load in one step, checks whether already compiled

# -------------------------------------------#
# Run TMB models
# -------------------------------------------#
mod_out <- lapply(X=models, FUN=run_model) # make a list of model outputs, run models
names(mod_out) <- models

mod_out # Aggregate_LRPs model has NaN for all Std. Errors

# Sgen values from Aggregate_LRPs model are much higher (at least double) those from 
# ricker_SMSY_Sgen. Need to look into why. 

# Try fixing values, optimize ricker parameters first

# -------------------------------------------#
# Save plots
# -------------------------------------------#

# Save plot of recruits~spawners by CU
plot_SR_by_CU(dat)

# Plot results of models
plot_compare_mods(mod1 = models[3], mod2= models[4])

# -------------------------------------------#
# Save model output
# -------------------------------------------#

# make a data frame with model estimates and CU names (works with ricker_SMSY_Sgen results)
CUcols <- rep(unique(dat$CU), 5) # make a vector of the CU names to bind with the estimates
resdf <- data.frame(parameter = names(res$value), value = res$value, CU= CUcols) # bind estimates with CU names
resdfw <- pivot_wider(resdf, names_from= parameter, values_from=value) # long to wide format
resdfw$alpha <- exp(resdfw$logA) # get alpha
resdfw$SMSY_80 <- 0.8 * resdfw$SMSY # get 80% of SMSY

res_sum <- resdfw %>% select(CU, alpha, Sgen, SMSY_80) %>% pivot_longer(cols=c(alpha, Sgen, SMSY_80)) # make summary table with same layout as report 
#write.csv(res_sum, "output/ricker_est_to_compare.csv") # write to csv
