# ----------------------------
# Master file to run analysis
# ----------------------------


# read in required libraries
library(TMB)

# source functions
source("functions.r")

# download data and save to .csv if not already present
source("download_data.r")

# read in data
dat <- read.csv("data_in/chum_SR_data.csv", stringsAsFactors = FALSE)

