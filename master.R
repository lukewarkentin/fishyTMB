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
models <- c("ricker_basic",   # basic Ricker model, 1 stock
               "ricker_multi_CUs", # Ricker with multiple stocks
               "ricker_SMSY_Sgen", # Ricker with SMSY and Sgen estimation (Hilborn and Walters)
               "Aggregate_LRPs") # Ricker, SMSY, Sgen and LRP binomial model of aggregate abundance over Sgen (1 phase optimization)

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
mod_out <- lapply(X=models, FUN=run_model, phases=1) # make a list of model outputs, run models
names(mod_out) <- models

mod_out # Aggregate_LRPs model has NaN for all Std. Errors. Did not coverge

# Sgen values from Aggregate_LRPs model are much higher (at least double) those from 
# ricker_SMSY_Sgen. Need to look into why. 

# -----------------#
# Try 2 phase optimization: fixing binomial parameters, optimize ricker parameters first
# -----------------#
map = list(logSgen = rep(factor(NA), length(unique(dat$CU))), B_0 = factor(NA), B_1 = factor(NA)) # fix logistic binomial model parameters
obj <- MakeADFun(data= model_input_list[["Aggregate_LRPs"]]$data_in, # make objective model function with fixed values of B_0 and B_1
                 parameters = model_input_list[["Aggregate_LRPs"]]$param_in, 
                 DLL="Aggregate_LRPs", map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr) # optimize
data.frame(summary(sdreport(obj)))
param1 <- obj$env$parList(opt$par) # get parameter estimates after phase 1 estimation

# pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj))) # get data frame of estimates with std error
All_Ests$Param <- row.names(All_Ests) # add column of parameter names
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ] # pull out vector of SMSY values
param1$logSgen <- log(0.3*SMSYs) # set initial value of logSgen as a function of SMSYs 

# phase 2 of optimization (for logSgen)
map2 <- list( B_0 = factor(NA), B_1 = factor(NA))
obj <- MakeADFun(data= model_input_list[["Aggregate_LRPs"]]$data_in, # make objective model function with fixed values of B_0 and B_1
                 parameters = param1,  # use parameter estimates from phase 1
                 DLL="Aggregate_LRPs", map=map2)

opt <- nlminb(obj$par, obj$fn, obj$gr) # optimize (phase 2)
data.frame(summary(sdreport(obj)))
param2 <- obj$env$parList(opt$par) # get parameter estimates after phase 2 estimation

# phase 3 of optimization (for B_0 and B_1)
obj <- MakeADFun(data= model_input_list[["Aggregate_LRPs"]]$data_in, # make objective model function with fixed values of B_0 and B_1
                 parameters = param2,  # use parameter estimates from phase 1
                 DLL="Aggregate_LRPs")

opt <- nlminb(obj$par, obj$fn, obj$gr) # optimize (phase 3) - try with control
data.frame(summary(sdreport(obj))) # B_0, B_1 and Agg_BM are not converging

mod_out$Aggregate_LRPs <- summary(sdreport(obj)) # replace 1 phase results with 3 phase. 

# -------------------------------------------#
# Save plots
# -------------------------------------------#

# Save plot of recruits~spawners by CU
plot_SR_by_CU(dat)

# Plot results of models, comparing ricker curves
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
