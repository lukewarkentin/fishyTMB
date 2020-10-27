# ----------------------------
# Master file to run stock recruit analysis
# ----------------------------

setwd("C:/github/fishyTMB") # temporary, set working directory

# read in required libraries
library(TMB)
library(ggplot2)

# source functions and figure functions
source("functions.r")
source("figures.r")

# download data and save to .csv if not already present, load and format data
source("load_data.r")

# Save plot of recruits~spawners by CU
plot_SR_by_CU(dat)





