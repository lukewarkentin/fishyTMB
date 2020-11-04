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

# -------------------------------------------#
# Prepare data for TMB models
# -------------------------------------------#
#model_input <- make_model_input(model_name = "ricker_basic", SRdat = dat)
#model_input <- make_model_input(model_name = "ricker_multi_CUs", SRdat=dat)
#model_input <- make_model_input(model_name = "ricker_SMSY_Sgen", SRdat=dat)
model_input <- make_model_input(model_name = "Aggregate_LRPs", SRdat=dat)

# -------------------------------------------#
# Compile TMB models
# -------------------------------------------#

# Compile and load TMB files
compile("TMB/ricker_basic.cpp") 
dyn.load(dynlib("TMB/ricker_basic"))

compile("TMB/ricker_multi_CUs.cpp") 
dyn.load(dynlib("TMB/ricker_multi_CUs"))

compile("TMB/ricker_SMSY_Sgen.cpp") 
dyn.load(dynlib("TMB/ricker_SMSY_Sgen"))

compile("TMB/Aggregate_LRPs.cpp")
dyn.load(dynlib("TMB/Aggregate_LRPs"))

# loadTMB("ricker_basic") # function to compile and load in one step, checks whether already compiled

# -------------------------------------------#
# Run TMB models
# -------------------------------------------#

# Make model function
#obj <- MakeADFun(data= model_input$data_in, parameters = model_input$param_in, DLL="ricker_basic")
#obj <- MakeADFun(data= model_input$data_in, parameters = model_input$param_in, DLL="ricker_multi_CUs")
#obj <- MakeADFun(data= model_input$data_in, parameters = model_input$param_in, DLL="ricker_SMSY_Sgen")
obj <- MakeADFun(data= model_input$data_in, parameters = model_input$param_in, DLL="Aggregate_LRPs")

# Optimize
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdreport(obj)
res <- sdreport(obj) # save results
round(res$value, 1)
# Sgen values from Aggregate_LRPs model are much higher (at least double) those from 
# ricker_SMSY_Sgen. Need to look into why.

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

# -------------------------------------------#
# Save plots
# -------------------------------------------#

# Save plot of recruits~spawners by CU
plot_SR_by_CU(dat)

# Plot results of ricker_SMSY_Sgen model
png("figures/fig_ricker_multi_CUs.png", height=400, width=700, units="px", pointsize=12 )
plot(exp(model_input$data_in$logR) ~ model_input$data_in$S, col = model_input$data_in$stock +1, xlab="Spawners", ylab="Recruits")
legend(x=700000, y = 1400000, legend=unique(dat$CU), col=1:7, pch=1)
# plot Ricker curves:
# R = alpha * S * exp(- beta * S)
for (i in 1:7) {
  curve(exp(res$value[grep("logA", names(res$value))][i]) * x * exp(-res$value[grep("B", names(res$value))][i] * x), 
        col=i, add=TRUE)
  # plot SMSY and Sgen
  abline(v=res$value[grep("SMSY", names(res$value))][i], col=i)
  abline(v=res$value[grep("Sgen", names(res$value))][i], col=i, lty=2)
}
dev.off()
