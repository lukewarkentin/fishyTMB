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
# 3 phase model gives convergence for B_0, B_1, or Agg_BM if data is scaled.

# Run 3 phase Aggregate LRP model with bounds on Sgen, and controls on optimization ---------

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
 upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than SMSY
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
 #upper[names(upper) =="B_0"] <- -2.3 # constrain B_0 to be less than -2.3
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
 CU_names <-  unique(dat$CU)
 mres$CU_name[!(mres$param %in% c("Agg_BM", "B_0", "B_1", "logit_preds"))] <- CU_names  # add a CU_name column. Reuses vector of CU names
 mres

 # # get predictions
 # preds <- inv_logit(mres$Estimate[mres$param=="logit_preds"])
 # preds_up <- inv_logit(mres$Estimate[mres$param=="logit_preds"] + mres$Std..Error[mres$param=="logit_preds"])
 # preds_low <- inv_logit(mres$Estimate[mres$param=="logit_preds"] - mres$Std..Error[mres$param=="logit_preds"])
 # # get the values to predict over
 # agg_abund <- model_input_list[["Aggregate_LRPs"]]$data_in$spawners_range * scale
 # 


# ------------------------#
# Get predictions for logistic model
# ------------------------#
preds <- inv_logit(mod_out$Aggregate_LRPs_2phase$Estimate[mod_out$Aggregate_LRPs_2phase$param=="logit_preds"])
preds_up <- inv_logit(mod_out$Aggregate_LRPs_2phase$Estimate[mod_out$Aggregate_LRPs_2phase$param=="logit_preds"] + mod_out$Aggregate_LRPs_2phase$Std..Error[mod_out$Aggregate_LRPs_2phase$param=="logit_preds"])
preds_low <- inv_logit(mod_out$Aggregate_LRPs_2phase$Estimate[mod_out$Aggregate_LRPs_2phase$param=="logit_preds"] - mod_out$Aggregate_LRPs_2phase$Std..Error[mod_out$Aggregate_LRPs_2phase$param=="logit_preds"])
# get the values to predict over
agg_abund <- model_input_list[["Aggregate_LRPs"]]$data_in$spawners_range * scale

# ----------------------#
# Plot logistic regression with benchmark
# ----------------------#
#png("figures/fig_check_logistic.png", width=8, height=6, units="in", pointsize=12, res=300)
# plot predicted values
plot( preds ~ agg_abund, type="l", ylim=c(0,1), lwd=2, col="dodgerblue", 
      xlab="Aggregate spawner abundance", ylab="proportion CUs > Sgen", 
      main="TMB model results, unconstrained B_0 (3 phase)")
# plot std error 
polygon(x=c(agg_abund, rev(agg_abund)), y= c(preds_up, rev(preds_low)), col=adjustcolor("dodgerblue", alpha=0.5), border=NA)

# plot benchmark
abline(v=mod_out$Aggregate_LRPs_2phase$Estimate[mod_out$Aggregate_LRPs_2phase$param=="Agg_BM"], col="gray", lty=2)
# # Get binomial regression parameters from above manual 3 phase with controls, etc.
# # plot observed data
points(y = obj$report()$N_Above_LRP/7, x= obj$report()$Agg_Abund*scale) # from manual 3 phase model fitting above
# B_0 <- mres$Estimate[mres$param=="B_0"]
# B_1 <- mres$Estimate[mres$param=="B_1"]
# plot binomial model estimate
# curve( inv_logit(B_0 + B_1*x/scale), col="red", lty=3 , lwd=2, add=TRUE)

# Get binomial regression parameters from other model outputs and plot
B_0 <- mod_out$Aggregate_LRPs$Estimate[mod_out$Aggregate_LRPs$param=="B_0"]
B_1 <- mod_out$Aggregate_LRPs$Estimate[mod_out$Aggregate_LRPs$param=="B_1"]
curve( inv_logit(B_0 + B_1*x/scale), col="orange", lty=3 , lwd=2, add=TRUE) # the one phase model results are different than the others. Does not fit data well
B_0 <- mod_out$Aggregate_LRPs_2phase$Estimate[mod_out$Aggregate_LRPs$param=="B_0"]
B_1 <- mod_out$Aggregate_LRPs_2phase$Estimate[mod_out$Aggregate_LRPs$param=="B_1"]
curve( inv_logit(B_0 + B_1*x/scale), col="pink", lty=3 , lwd=2, add=TRUE) # 2 phase model gives same results
B_0 <- mod_out$Aggregate_LRPs_3phase$Estimate[mod_out$Aggregate_LRPs$param=="B_0"]
B_1 <- mod_out$Aggregate_LRPs_3phase$Estimate[mod_out$Aggregate_LRPs$param=="B_1"]
curve( inv_logit(B_0 + B_1*x/scale), col="brown", lty=3 , lwd=2, add=TRUE) # 3 phase model gives same results 

#legend(x=1000000, y=0.2, legend=c("Predicted", "Formula", "benchmark, p=0.8"), lty=c(1,2,2), lwd=c(2,3,1),col=c("dodgerblue", "red", "gray"), )
#dev.off()

# --------------------------------#
# Use glm function to compare different models
# --------------------------------#

#### Compare with and without influential point, logit vs. cauchit link functions

# make a data frame of just aggregate abundance and prop over Sgen
tdf <- data.frame(prop = obj$report()$N_Above_LRP/7 , abd = obj$report()$Agg_Abund) 
tdf2 <- tdf[-which.max(tdf$abd),] # remove highest abundance, influential point

#### Compare logit to cauchit fits (keep influential point)
fit_l <- glm(prop ~ abd, data=tdf, family=quasibinomial(link="logit"))
fit_l2 <- glm(prop ~ abd, data=tdf2, family=quasibinomial(link="logit")) # logit link without influential point
fit_c <- glm(prop ~ abd, data=tdf, family=quasibinomial(link="cauchit"))
fit_c2 <- glm(prop ~ abd, data=tdf2, family=quasibinomial(link="cauchit")) # cauchit link without influential point
# fit_p <- glm(prop ~ abd, data=tdf, family=quasibinomial(link="probit")) 
# fit_p2 <- glm(prop ~ abd, data=tdf2, family=quasibinomial(link="probit"))
# Note: probit link gives an intermediate curve under inflection point, 
# and higher than either logit or cauchit above inflection point.

coefs_l <- fit_l$coefficients
coefs_l2 <- fit_l2$coefficients
coefs_c <- fit_c$coefficients
coefs_c2 <- fit_c2$coefficients
# coefs_p <- fit_p$coefficients
# coefs_p2 <- fit_p2$coefficients

png("figures/fig_logit_vs_cauchit.png", width=8, height=6, units="in", pointsize=12, res=300)
plot(tdf$prop ~ tdf$abd, xlim=c(-100, 300), ylim=c(0,1), xlab="aggregate abundance (scaled)", ylab="proportion CU>Sgen", 
     main="Generic glm() results")
points(tdf$prop[50] ~ tdf$abd[50], xlim=c(-100, 300), ylim=c(0,1), col="red")
curve(inv_logit(coefs_l[1] + coefs_l[2] * x ), add=TRUE)
curve(inv_logit(coefs_l2[1] + coefs_l2[2] * x ), add=TRUE, col="gray")
curve(VGAM::cauchitlink(coefs_c[1] + coefs_c[2] * x , inverse=TRUE), add=TRUE, col="purple") # inverse cauchit link from VGAM package
curve(VGAM::cauchitlink(coefs_c2[1] + coefs_c2[2] * x , inverse=TRUE), add=TRUE, col="pink") # inverse cauchit link from VGAM package
# curve(VGAM::probitlink(coefs_p[1] + coefs_p[2] * x, inverse=TRUE), add=TRUE, col="orange") # inverse probit link from VGAM package
# curve(VGAM::probitlink(coefs_p2[1] + coefs_p2[2] * x, inverse=TRUE), add=TRUE, col="orange3") # inverse probit link from VGAM package
legend(x=100, y=0.4, lty=1, col=c("black", "gray", "purple", "pink"), legend=c("logit", "logit (no influential point)", "cauchit", "cauchit (no influential point)"))
dev.off()

# -------------------------------------------#
# Save plots
# -------------------------------------------#

# Plot results of models, comparing ricker curves # need to unscale parameter outputs for this to work
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

# plot spread in spawner abundances
# compare to Coho data
options(scipen=1000000)

coho_dat <- read.csv("https://raw.githubusercontent.com/Pacific-salmon-assess/SalmonLRP_RetroEval/master/IFCohoStudy/DataIn/IFCoho_SRbyCU.csv") # read in coho data from local repository
# rename columns to match chum data
names(coho_dat)[grep("CU_Name", names(coho_dat))] <- "CU"
names(coho_dat)[grep("Spawners", names(coho_dat))] <- "spawners"
names(coho_dat)[grep("Recruits", names(coho_dat))] <- "recruits"
names(coho_dat)[grep("BroodYear", names(coho_dat))] <- "year"
# select data needed
coho_dat <- coho_dat[,c(1,3,4,5)]
coho_dat$sp <- "coho" # add species column
dat$sp <- "chum"
all_dat <- rbind(dat, coho_dat) # bind rows
# plot comparison of distributions
png(filename="figures/fig_compare_chum_coho_spawner_dist.png", width=8, height=6,units="in", res=300)
ggplot(all_dat, aes(x=spawners, colour=CU, fill=CU)) +
  geom_point(aes(y=0, x=spawners), shape=108, colour="black", size=2) +
  geom_density(alpha=0.5) + 
  scale_x_log10( breaks= c(10^(1:10))) +
  facet_wrap(~sp, ncol=1) +
  xlab("log10 spawners") +
  ylab("Density") +
  theme_classic()
dev.off()

# compare correlation among CUs
# convert long to wide data for correlations
coho_dat_w <- coho_dat %>% select(CU, spawners, year) %>% pivot_wider(names_from=CU, values_from=spawners)
chum_dat_w <- dat %>% select(CU, spawners, year) %>% pivot_wider(names_from=CU, values_from=spawners)
# Plot correlation of spawner abundances
png(filename="figures/fig_cor_spawners_coho.png", width=6, height=6, units="in", res=300)
PerformanceAnalytics::chart.Correlation(coho_dat_w[,-1])
dev.off()

png(filename="figures/fig_cor_spawners_chum.png", width=6, height=6, units="in", res=300)
PerformanceAnalytics::chart.Correlation(chum_dat_w[,-1])
dev.off()


# -------------------------------------------#
# Save model output (to compare with Holt et al. 2018)
# -------------------------------------------#

# make a data frame with model estimates and CU names (works with any of the mod_out models)
resdf <- mod_out$Aggregate_LRPs_2phase # bind estimates with CU names
resdfw <- pivot_wider(resdf[!is.na(resdf$CU_name) , -grep("Std..Error", names(resdf))], names_from = param, values_from=Estimate) # long to wide format
resdfw$SMSY_80 <- 0.8 * resdfw$SMSY # get 80% of SMSY

res_sum <- resdfw %>% select(CU_name, A, Sgen, SMSY_80) %>% pivot_longer(cols=c(A, Sgen, SMSY_80)) # make summary table with same layout as report 
#write.csv(res_sum, "output/ricker_est_to_compare.csv") # write to csv

# Look at different distributions
all_dat_y <- all_dat %>% group_by(year, sp) %>% summarise(agg_abd = sum(spawners))
hist(all_dat_y[all_dat_y$sp=="coho",]$agg_abd, breaks=10)

x<- seq(0, max(all_dat_y[all_dat_y$sp=="chum",]$agg_abd*1.1), 100)
yc <- dcauchy(x, location=mean(all_dat_y[all_dat_y$sp=="chum",]$agg_abd), scale=100000)
yn <- dnorm(x, mean=mean(all_dat_y[all_dat_y$sp=="chum",]$agg_abd), sd=100000)
yl <- dlogis(x, location= mean(all_dat_y[all_dat_y$sp=="chum",]$agg_abd), scale=100000)
png("figures/fig_compare_distributions.png", width=6, height=5, units="in", res=300)
plot(x=x, y=yc, type="l", col="dodgerblue", xlim=c(0, max(all_dat_y[all_dat_y$sp=="chum",]$agg_abd)*1.05), 
     xlab="Aggregage spawner abundance", ylab="density", lwd=2)
hist(all_dat_y[all_dat_y$sp=="chum",]$agg_abd, freq=FALSE, breaks=10, add=TRUE, col=NULL)
lines(x=x, y=yn, type="l", col="purple", lwd=2)
lines(x=x, y=yl, type="l", col="orange", lwd=2)
legend(x=1500000, y=0.000002, col=c("dodgerblue", "purple", "orange"), legend=c("Cauchy", "Normal", "Logistic"),lwd=2, lty=1)
dev.off()

