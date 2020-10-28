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
compile("TMB/ricker.cpp","-O1 -g",DLLFLAGS="") # extra arguments allow debug with gdbsource()
dyn.load(dynlib("TMB/ricker"))
# loadTMB("ricker") # function to compile and load in one step, checks whether already compiled

# TMB input data (must be list)
data <- list() # create empty list
data$S <- dat$spawners # spawners 
data$logR <- log(dat$recruits) # natural log of recruits
data$stock <- as.integer(as.factor(dat$CU)) # numeric vector of CU (conservation unit)
n_stocks <- length(unique(dat$CU)) # number of CUs
data$n_stocks <- n_stocks

# Set-up parameter list for TMB function calls
param <- list()
param$logA <- rep(1, n_stocks) # intial values of logA for each CU
param$logB <- as.numeric(log(1/( (dat %>% group_by(CU) %>% summarise(x=quantile(spawners, 0.8)))$x) ))
param$logSigma <- rep(-2, n_stocks)

# Make model function
obj <- MakeADFun(data, param, DLL="ricker") # This makes R abort and give this error: 
# TMB has received an error from EIgen. The following conditions was not met:
# index >=0 && index < size()
# Please chekc your matrix-vector bounds etc., or run your program through a debugger
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt
sdreport(obj)

# check for out of bounds error with debug function gdbsource()
gdbsource("ricker.R") # gives: 
#Error in system(cmd, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE) : 
#'gdb' not found