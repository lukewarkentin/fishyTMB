# -----------------------------------
# Functions
#
#
#   get_SR_dat : Test if data is already present; if false, write data to csv
#   loadTMB : Test whether TMB files are already compiled, if not, don't recompile. Then load model
#   make_model_input: Prepare data and parameter inputs for a particular TMB model
#
# -----------------------------------

# Function to test if data is already present; if false, write data to csv
get_SR_dat <- function(x) {
  if(file.exists(paste0("data_in/", x))) {  # test if data is already present in folder
          warning("File already exists") # if file exists, give warning message
         # if file does not exist, download .csv file from SalmonLRP_RetroEval github repository
      } else { dat <- read.csv("https://raw.githubusercontent.com/Pacific-salmon-assess/SalmonLRP_RetroEval/withChum/SCChumStudy/DataOut/SRdatWild.csv", stringsAsFactors = FALSE)
         dat1 <- dat[, names(dat) %in% c("Year", "CU", "Escape", "Recruit")] # select columns to keep
         names(dat1) <- c("year", "CU", "spawners", "recruits") # rename columns
         write.csv(dat1, paste0("data_in/", x), row.names = FALSE)
      } # end of else statement
}

# Function to test whether TMB files are already compiled, if not, don't recompile. Then load model
loadTMB <- function(x) {
  modpath <- paste0("TMB/", x) # model path string
  if(file.exists(paste0(modpath, ".dll"))) { # test if TMB .cpp file is already present
        warning("Already compiled: .dll file already exists") # if file exists, give warning
      } else { compile(paste0(modpath, ".cpp"), "-O1 -g",DLLFLAGS="") # compile .cpp file of model, extra arguments allow debug with gdbsource()
  } # end of else statement 
  dyn.load(dynlib(modpath)) # load .dll file of model
}

# Function to prepare data and parameter inputs for a particular TMB model
make_model_input <- function(model_name, SRdat) {
  data_in <- list() # create empty list for observed data to go into model
  param_in <- list() # create empty list for parameters with initial values assigned
  
  # if model is basic Ricker with 1 population
  if(model_name =="ricker_basic") {
    SRdat <- SRdat[grep("1 -", SRdat$CU), ] # if model is the basic ricker with one CU, grab just one CU
    # data
    data_in$S <- SRdat$spawners # spawners 
    data_in$logR <- log(SRdat$recruits) # natural log of recruits 
    # parameters
    param_in$logA <- 1
    param_in$logB <- as.numeric(log(quantile(SRdat$spawners, 0.8)))
    param_in$logSigma <- -2
  }
  # Set spawner and log(recruit) data inputs

    
  # If model is not basic ricker (multiple populations):
  if(model_name=="ricker_multi_CUs") {
    # data
    data_in$S <- SRdat$spawners # spawners 
    data_in$logR <- log(SRdat$recruits) # natural log of recruitsparam_in$logA <- 1
    data_in$stock <- as.integer(as.factor(SRdat$CU)) # numeric vector of CU (conservation unit)
    n_stocks <- length(unique(SRdat$CU)) # number of CUs
    #data_in$n_stocks <- n_stocks # number of stocks
    # parameters
    param_in$logA <- rep(1, n_stocks) # intial values of logA for each CU
    param_in$logB <- as.numeric(log(1/( (SRdat %>% group_by(CU) %>% summarise(x=quantile(spawners, 0.8)))$x) ))
    param_in$logSigma <- rep(-2, n_stocks)
  }
  
  model_in <- list()
  model_in$data_in <- data_in
  model_in$param_in <- param_in
  
  model_in
  
} # end of function

