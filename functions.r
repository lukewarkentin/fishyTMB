# -----------------------------------
# Functions
#
#   inv_logit
#   logit
#   get_SR_dat : Test if data is already present; if false, write data to csv
#   make_model_input: Prepare data and parameter inputs for a particular TMB model and data
#   run_model: Run TMB model and get output
#
# -----------------------------------

# Inverse logit and logit funcitons can come in handy =====================================================================
inv_logit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  log(x/1-x)
}


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

# Function to prepare data and parameter inputs for a particular TMB model
make_model_input <- function(model_name, SRdat) {
  
  data_in <- list() # create empty list for observed data to go into model
  param_in <- list() # create empty list for parameters with initial values assigned
  
  # if model is basic Ricker with 1 population
  if(model_name =="ricker_basic") {
    SRdat <- SRdat[grep("1 -", SRdat$CU), ] # if model is the basic ricker with one CU, grab just one CU
    # parameters (only one CU)
    param_in$logA <- 1
    param_in$logB <- as.numeric(log(quantile(SRdat$spawners, 0.8)))
    param_in$logSigma <- -2
  }
  
  # Set spawner and log(recruit) data inputs (for all models)
    data_in$S <- SRdat$spawners # spawners 
    data_in$logR <- log(SRdat$recruits) # natural log of recruits 
    
  # If multi-CU model (all models except basic)
  if(model_name != "ricker_basic") {
    # data
    data_in$stock <- as.integer(as.factor(SRdat$CU)) - 1 # numeric vector of CU (conservation unit)
    n_stocks <- length(unique(SRdat$CU)) # number of CUs
    data_in$n_stocks <- n_stocks # number of stocks
    # parameters
    param_in$logA <- rep(1, n_stocks) # intial values of logA for each CU
    param_in$logB <- as.numeric(log(1/( (SRdat %>% group_by(CU) %>% summarise(x=quantile(spawners, 0.8)))$x) ))
    param_in$logSigma <- rep(-2, n_stocks)
  }
  
  # If model is with SMSY and Sgen outputs (all models except ricker_basic and ricker_multi_CUs)
  if(! model_name %in% c("ricker_basic", "ricker_multi_CUs")) {
    # data
    # parameters
    param_in$logSgen <- log((SRdat %>% group_by(CU) %>%  summarise(x=quantile(spawners, 0.5)))$x) 
  }
  
  # If model is has logistic regression
  if(model_name == "Aggregate_LRPs") {
    # data
    data_in$yr <- SRdat$year
    data_in$Mod_Yr_0 <- min(SRdat$year)
    data_in$Mod_Yr_n <- max(SRdat$year)
    agg_data <- SRdat %>% group_by(year) %>% summarise(total_spawners = sum(spawners, na.rm=TRUE))
    spawners_range <- seq(0,max(agg_data$total_spawners),length.out = 100)
    data_in$spawners_range <- spawners_range # vector to predict N over threshold
    # parameters
    param_in$B_0 <- 2 # from Brooke's older code
    param_in$B_1 <- 0.1 # from Brooke's older code
  }
  
  model_in <- list()
  model_in$data_in <- data_in
  model_in$param_in <- param_in
  
  model_in
  
} # end of function


# Function to run optimization and save output
run_model <- function(model_name, phases, CU_names) {
  
  # single phase model =================
  
  if(phases == 1 ){ # for single phase model
    obj <- MakeADFun(data= model_input_list[[model_name]]$data_in, parameters = model_input_list[[model_name]]$param_in, DLL=model_name)
    # Optimize
    opt <- nlminb(obj$par, obj$fn, obj$gr)
  } 
  
  # 2 phase model =================== 
  # 1 phase and 3 phase ricker estimates were different. See if optimization Sgen with alpha and beta makes a difference
  
  if(phases==2) {
    map = list(B_0 = factor(NA), B_1 = factor(NA)) # fix logistic binomial model parameters
    obj <- MakeADFun(data= model_input_list[[model_name]]$data_in, # make objective model function with fixed values of B_0 and B_1
                     parameters = model_input_list[[model_name]]$param_in, 
                     DLL=model_name, map=map)
    opt <- nlminb(obj$par, obj$fn, obj$gr) # optimize
    param1 <- obj$env$parList(opt$par) # get parameter estimates after phase 1 estimation
    
    # phase 2 of optimization (B_0 and B_1)
    obj <- MakeADFun(data= model_input_list[[model_name]]$data_in,
                     parameters = param1,  # use parameter estimates from phase 1
                     DLL=model_name)
    
    opt <- nlminb(obj$par, obj$fn, obj$gr) # optimize (phase 2)
  }
  
  # 3 phase model ===================
  
  if(phases == 3) {
    map = list(logSgen = rep(factor(NA), length(unique(dat$CU))), B_0 = factor(NA), B_1 = factor(NA)) # fix logistic binomial model parameters
    obj <- MakeADFun(data= model_input_list[[model_name]]$data_in, # make objective model function with fixed values of B_0 and B_1
                     parameters = model_input_list[[model_name]]$param_in, 
                     DLL=model_name, map=map)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5)) # optimize, try control
    param1 <- obj$env$parList(opt$par) # get parameter estimates after phase 1 estimation
    
    # pull out SMSY values
    All_Ests <- data.frame(summary(sdreport(obj))) # get data frame of estimates with std error
    All_Ests$Param <- row.names(All_Ests) # add column of parameter names
    SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ] # pull out vector of SMSY values
    param1$logSgen <- log(0.3*SMSYs) # set initial value of logSgen as a function of SMSYs 
    
    # phase 2 of optimization (for logSgen)
    map2 <- list( B_0 = factor(NA), B_1 = factor(NA))
    obj <- MakeADFun(data= model_input_list[[model_name]]$data_in, # make objective model function with fixed values of B_0 and B_1
                     parameters = param1,  # use parameter estimates from phase 1
                     DLL=model_name, map=map2)
    
    opt <- nlminb(obj$par, obj$fn, obj$gr) # optimize (phase 2)
    param2 <- obj$env$parList(opt$par) # get parameter estimates after phase 2 estimation
    
    # phase 3 of optimization (for B_0 and B_1)
    obj <- MakeADFun(data= model_input_list[[model_name]]$data_in, # make objective model function with fixed values of B_0 and B_1
                     parameters = param2,  # use parameter estimates from phase 1
                     DLL=model_name)
    opt <- nlminb(obj$par, obj$fn, obj$gr) # optimize (phase 3)
    }
  
    # make data frame of model results ================
  
    mres <- data.frame(summary(sdreport(obj))) 
    mres$param <- row.names(mres) # make column of parameter names
    mres$param <- sub("\\.\\d*", "", mres$param ) # remove .1, .2 etc from parameter names. \\. is "." and \\d means any digit
    mres$mod <- model_name # this would be to add a column of model names. Would work if in a function
    mres$phases <- phases # number of phases
    mres$CU_ID[!(mres$param %in% c("Agg_BM", "B_0", "B_1", "logit_preds"))] <- seq_along(CU_names)  # add a CU_ID column
    mres <- merge(mres, data.frame(CU_name = CU_names, CU_ID= seq_along(CU_names)), by="CU_ID", all.x=TRUE) # merge CU names
    mres <- mres[order(mres$param),] # order based on parameter
    mres
    
} # end of function
