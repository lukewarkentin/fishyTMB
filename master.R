# ----------------------------
# Master file to run stock recruit analysis
# ----------------------------

setwd("C:/github/fishyTMB") # temporary, set working directory

# read in required libraries
library(TMB)
library(ggplot2)
library(dplyr)

# source functions and figure functions
source("functions.r")
source("figures.r")

# download data and save to .csv if not already present, load and format data
source("load_data.r")

# Save plot of recruits~spawners by CU
plot_SR_by_CU(dat)

# Compile and load TMB files if not done already

compile("TMB/ricker_basic.cpp","-O1 -g",DLLFLAGS="") # extra arguments allow debug with gdbsource()
dyn.load(dynlib("TMB/ricker_basic"))

compile("TMB/ricker_multi_CUs.cpp") # extra arguments allow debug with gdbsource()
dyn.load(dynlib("TMB/ricker_multi_CUs"))

# loadTMB("ricker_basic") # function to compile and load in one step, checks whether already compiled

# Prepare model input data
#model_input <- make_model_input(model_name = "ricker_basic", SRdat = dat)
model_input <- make_model_input(model_name = "ricker_multi_CUs", SRdat=dat)

# Make model function
#obj <- MakeADFun(data= model_input$data_in, parameters = model_input$param_in, DLL="ricker_basic")
obj <- MakeADFun(data= model_input$data_in, parameters = model_input$param_in, DLL="ricker_multi_CUs")

# Optimize
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdreport(obj)
res <- sdreport(obj)
res$value
# plot results to check
plot(data$logR ~ data$S)
# plot formula:
# logR = logA + logS - BS
curve(res$value["logA"] + log(x) - res$value["B"] * x , add=TRUE)

# check for out of bounds error with debug function gdbsource()
gdbsource("ricker.R") # gives: 
#Error in system(cmd, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE) : 
#'gdb' not found