# -----------------------------------
# Read in data and save 
# -----------------------------------

# Download South Coast Chum stock-recruit data
# Raw data provided by Pieter Van Will, infilling code by Brooke Davis (see https://github.com/Pacific-salmon-assess/SalmonLRP_RetroEval.git)
# save data to csv if not already present
get_SR_dat("chum_SR_data.csv")

# read in data
dat <- read.csv("data_in/chum_SR_data.csv")
dat <- dat[!is.na(dat$recruits), ] # remove NA rows