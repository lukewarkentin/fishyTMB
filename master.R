# Master file to run stock recruit analysis

# -------------------------------------------#
# Load packages, data, and functions
# -------------------------------------------#

# setwd("C:/github/fishyTMB") # temporary, set working directory

# read in required libraries
library(TMB)
library(ggplot2)
library(dplyr)

# source functions and figure functions
source("functions.r")
source("figures.r")

# download data and save to .csv if not already present, load and format data
source("load_data.r")

# -------------------------------------------#
# Prepare data for TMB models
# -------------------------------------------#
#model_input <- make_model_input(model_name = "ricker_basic", SRdat = dat)
model_input <- make_model_input(model_name = "ricker_multi_CUs", SRdat=dat)

# -------------------------------------------#
# Compile TMB models
# -------------------------------------------#

# Compile and load TMB files
compile("TMB/ricker_basic.cpp") 
dyn.load(dynlib("TMB/ricker_basic"))

compile("TMB/ricker_multi_CUs.cpp") 
dyn.load(dynlib("TMB/ricker_multi_CUs"))

# loadTMB("ricker_basic") # function to compile and load in one step, checks whether already compiled

# -------------------------------------------#
# Run TMB models
# -------------------------------------------#

# Make model function
#obj <- MakeADFun(data= model_input$data_in, parameters = model_input$param_in, DLL="ricker_basic")
obj <- MakeADFun(data= model_input$data_in, parameters = model_input$param_in, DLL="ricker_multi_CUs")

# Optimize
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdreport(obj)
res <- sdreport(obj)
res$value

# -------------------------------------------#
# Save plots
# -------------------------------------------#

# Save plot of recruits~spawners by CU
plot_SR_by_CU(dat)

# plot results to check
png("figures/fig_ricker_multi_CUs.png", height=400, width=700, units="px", pointsize=12 )
plot(exp(model_input$data_in$logR) ~ model_input$data_in$S, col = model_input$data_in$stock +1, xlab="Spawners", ylab="Recruits")
legend(x=700000, y = 1400000, legend=unique(dat$CU), col=1:7, pch=1)
# plot Ricker curves:
# R = alpha * S * exp(- beta * S)
for (i in 1:7) {
  curve(exp(res$value[grep("logA", names(res$value))][i]) * x * exp(-res$value[grep("B", names(res$value))][i] * x), 
        col=i, add=TRUE)
}
dev.off()
