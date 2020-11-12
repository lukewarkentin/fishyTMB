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

scale<-10000 # set scale at 10000. TMB doesn't work as well with very large numbers

# make list of model input lists
model_input_list <- lapply(X=models, make_model_input, SRdat = dat, scale=scale) # model_name argument automatically evaluates to elements of models object
names(model_input_list) <- models # name elements of list after model names

# -------------------------------------------#
# Run TMB models
# -------------------------------------------#
mod_out <- lapply(X=models, FUN=run_model, phases=1, CU_names = unique(dat$CU)) # make a list of model outputs, run models
names(mod_out) <- models
mod_out # Aggregate_LRPs model has NaN for all Std. Errors. Did not coverge

# Run 2 and 3 phase models, append to model results list
mod_out$Aggregate_LRPs_2phase <- run_model(model_name = "Aggregate_LRPs", phases=2, CU_names=unique(dat$CU))
mod_out$Aggregate_LRPs_3phase <- run_model(model_name = "Aggregate_LRPs", phases=3, CU_names=unique(dat$CU))
mod_out
# 3 phase model gives convergence for ricker parameters, but not B_0, B_1, or Agg_BM

# =======================================================================#
# --- Figure out how logistic predictions are saved from model object ---
# =======================================================================#
# Use Aggregate_LRP model
# 3 phase optimization

# Phase 1:

map = list(logSgen=factor(rep(NA, 7)),B_0 = factor(NA), B_1 = factor(NA)) # fix logSgen and logistic parameters

obj <- MakeADFun(data= model_input_list[["Aggregate_LRPs"]]$data_in, # make objective model function with fixed values of B_0 and B_1
                   parameters = model_input_list[["Aggregate_LRPs"]]$param_in, 
                   DLL="Aggregate_LRPs", map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr) # optimize
param1 <- obj$env$parList(opt$par) # get parameter estimates after phase 1 estimation




# pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
# set initial Sgen param as a function of Smsy
param1$logSgen <- log(0.3*SMSYs)


# Phase 2:
map = list(B_0 = factor(NA), B_1 = factor(NA)) # fix logistic binomial model parameters
obj <- MakeADFun(data= model_input_list[["Aggregate_LRPs"]]$data_in, # make objective model function with fixed values of B_0 and B_1
                 parameters = param1, 
                 DLL="Aggregate_LRPs", map=map)

## Create upper & lower bounds vectors that are same length and order as nlminb start vector
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.001) # constrain Sgen to be positive
lower<-unname(lower)


opt <- nlminb(obj$par, obj$fn, obj$gr,control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper, lower=lower) # optimize
param2 <- obj$env$parList(opt$par) # get parameter estimates after phase 1 estimation

  
# phase 3 of optimization (B_0 and B_1)
obj <- MakeADFun(data= model_input_list[["Aggregate_LRPs"]]$data_in,
                   parameters = param2,  # use parameter estimates from phase 2
                   DLL="Aggregate_LRPs")
  

## Create upper & lower bounds vectors that are same length and order as nlminb start vector
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy
upper[names(upper) =="B_0"] <- -2.3 # constrain B_0 to be less than -2.3
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<- -Inf
lower[names(lower) =="logSgen"] <- log(0.001) # constrain Sgen to be positive
lower<-unname(lower)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper, lower=lower) # optimize (phase 3)

# make data frame of results, and format
mres <- data.frame(summary(sdreport(obj))) 
mres$param <- row.names(mres) # make column of parameter names
mres$param <- sub("\\.\\d*", "", mres$param ) # remove .1, .2 etc from parameter names. \\. is "." and \\d means any digit

# This section messes up order of predictions
#CU_names <-  unique(dat$CU)
#mres$CU_ID[!(mres$param %in% c("Agg_BM", "B_0", "B_1", "logit_preds"))] <- seq_along(CU_names)  # add a CU_ID column
#mres <- merge(mres, data.frame(CU_name = CU_names, CU_ID= seq_along(CU_names)), by="CU_ID", all.x=TRUE) # merge CU names
mres


# get predictions 
preds <- inv_logit(mres$Estimate[mres$param=="logit_preds"])
preds_up <- inv_logit(mres$Estimate[mres$param=="logit_preds"] + mres$Std..Error[mres$param=="logit_preds"])
preds_low <- inv_logit(mres$Estimate[mres$param=="logit_preds"] - mres$Std..Error[mres$param=="logit_preds"])
# get the values to predict over
agg_abund <- model_input_list[["Aggregate_LRPs"]]$data_in$spawners_range * scale

# ----------------------#
# Plot logistic regression with benchmark
# ----------------------#
png("figures/fig_check_logistic.png", width=8, height=6, units="in", pointsize=12, res=300)
# plot predicted values
plot( preds ~ agg_abund, type="l", ylim=c(0,1), lwd=2, col="dodgerblue", xlab="Aggregate spawner abundance", ylab="proportion CUs > Sgen")
# plot std error 
polygon(x=c(agg_abund, rev(agg_abund)), y= c(preds_up, rev(preds_low)), col=adjustcolor("dodgerblue", alpha=0.5), border=NA)
# plot observed data
points(y = obj$report()$N_Above_LRP/7, x= obj$report()$Agg_Abund*scale)
# plot benchmark
abline(v=mres$Estimate[mres$param=="Agg_BM"]*scale, col="gray", lty=2)
# Get binomial regression parameters
B_0 <- mres$Estimate[mres$param=="B_0"]
B_1 <- mres$Estimate[mres$param=="B_1"]
# plot binomial model estimate
curve( inv_logit(B_0 + B_1*x/scale), col="red", lty=3 , lwd=2, add=TRUE)
legend(x=1000000, y=0.2, legend=c("Predicted", "Formula", "benchmark, p=0.8"), lty=c(1,2,2), lwd=c(2,3,1),col=c("dodgerblue", "red", "gray"), )
dev.off()

# -------------------------------------------#
# Save plots
# -------------------------------------------#

# Plot results of models, comparing ricker curves
plot_compare_mods(mod1 = models[2], mod2= models[3])
plot_compare_mods(mod1 = models[2], mod2= "Aggregate_LRPs_2phase")
plot_compare_mods(mod1 = models[3], mod2= "Aggregate_LRPs_2phase")
plot_compare_mods(mod1 = models[2], mod2= "Aggregate_LRPs_3phase")
plot_compare_mods(mod1 = models[3], mod2= "Aggregate_LRPs_3phase")


# Plot recruits/spawner over time
ggplot(dat, aes(y=recruits/spawners, x=year, colour=CU)) + 
  geom_line() + 
  geom_point() +
  scale_y_log10() +
  geom_hline(aes(yintercept=1)) +
  facet_wrap(~CU, scales="free_y")  
  #theme_classic()

ggplot(dat[dat$CU=="3 - Upper Knight", ], aes(y=recruits, x=spawners)) + 
  geom_point(size=4,aes(colour=year)) + 
  geom_text(aes(label=year)) + 
  theme_classic()

# -------------------------------------------#
# Save model output (to compare with Holt et al. 2018)
# -------------------------------------------#

# make a data frame with model estimates and CU names (works with ricker_SMSY_Sgen results)
CUcols <- rep(unique(dat$CU), 5) # make a vector of the CU names to bind with the estimates
resdf <- data.frame(parameter = names(res$value), value = res$value, CU= CUcols) # bind estimates with CU names
resdfw <- pivot_wider(resdf, names_from= parameter, values_from=value) # long to wide format
resdfw$alpha <- exp(resdfw$logA) # get alpha
resdfw$SMSY_80 <- 0.8 * resdfw$SMSY # get 80% of SMSY

res_sum <- resdfw %>% select(CU, alpha, Sgen, SMSY_80) %>% pivot_longer(cols=c(alpha, Sgen, SMSY_80)) # make summary table with same layout as report 
#write.csv(res_sum, "output/ricker_est_to_compare.csv") # write to csv
