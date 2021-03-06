# Master file to run stock recruit analysis 

# -------------------------------------------#
# Load packages, data, and functions
# -------------------------------------------#

setwd("C:/github/fishyTMB") # temporary, set working directory

# read in required libraries
library(TMB)
library(ggplot2)
options(scipen=1000000)
library(dplyr)
library(tidyr)
library(corrplot)

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
# Run TMB models with automated function
# -------------------------------------------#
mod_out <- lapply(X=models, FUN=run_model, phases=1, CU_names = unique(dat$CU)) # make a list of model outputs, run models
names(mod_out) <- models
mod_out # Aggregate_LRPs model has NaN for all Std. Errors. Did not coverge

# Run 2 and 3 phase models, append to model results list
mod_out$Aggregate_LRPs_2phase <- run_model(model_name = "Aggregate_LRPs", phases=2, CU_names=unique(dat$CU))
mod_out$Aggregate_LRPs_2phase <- run_model(model_name = "Aggregate_LRPs", phases=2, CU_names=unique(dat$CU))

mod_out$Aggregate_LRPs_3phase <- run_model(model_name = "Aggregate_LRPs", phases=3, CU_names=unique(dat$CU))
mod_out
# 3 phase model gives convergence for B_0, B_1, or Agg_BM if data is scaled.

# ----------------#
# Run 3 phase Aggregate LRP model with bounds on Sgen, and controls on optimization ---------
# ----------------#

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
# upper<-unlist(obj$par)
#  upper[1:length(upper)]<-Inf
#  upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than SMSY
#  upper<-unname(upper)
# 
#  lower<-unlist(obj$par)
#  lower[1:length(lower)]<--Inf
#  lower[names(lower) =="logSgen"] <- log(0.001) # constrain Sgen to be positive
#  lower<-unname(lower)


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
 #upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy
 upper[names(upper) =="B_0"] <- -2.3 # constrain B_0 to be less than -2.3
 upper<-unname(upper)

 lower<-unlist(obj$par)
 lower[1:length(lower)]<- -Inf
 #lower[names(lower) =="logSgen"] <- log(0.001) # constrain Sgen to be positive
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

#-------------
# ------------------------#
# Get predictions for logistic model
# ------------------------#
# From unconstrained B_0
preds <- inv_logit(mod_out$Aggregate_LRPs_3phase$Estimate[mod_out$Aggregate_LRPs_3phase$param=="logit_preds"])
preds_up <- inv_logit(mod_out$Aggregate_LRPs_3phase$Estimate[mod_out$Aggregate_LRPs_3phase$param=="logit_preds"] + mod_out$Aggregate_LRPs_3phase$Std..Error[mod_out$Aggregate_LRPs_3phase$param=="logit_preds"])
preds_low <- inv_logit(mod_out$Aggregate_LRPs_3phase$Estimate[mod_out$Aggregate_LRPs_3phase$param=="logit_preds"] - mod_out$Aggregate_LRPs_3phase$Std..Error[mod_out$Aggregate_LRPs_3phase$param=="logit_preds"])


# From B_0 constrained 
preds_c <- inv_logit(mres$Estimate[mres$param=="logit_preds"])
preds_c_up <- inv_logit(mres$Estimate[mres$param=="logit_preds"] + mres$Std..Error[mres$param=="logit_preds"])
preds_c_low <- inv_logit(mres$Estimate[mres$param=="logit_preds"] - mres$Std..Error[mres$param=="logit_preds"])

# get the values to predict over
agg_abund <- model_input_list[["Aggregate_LRPs"]]$data_in$spawners_range * scale

# ----------------------#
# Plot logistic regression with benchmark
# ----------------------#
png("figures/fig_check_logistic_base.png", width=8, height=6, units="in", pointsize=12, res=300)
# plot predicted values
plot( preds ~ agg_abund, type="l", ylim=c(0,1), xlim=c(0, max(agg_abund)*1.01), lwd=2, col="dodgerblue", 
      xlab="Aggregate spawner abundance", ylab="proportion CUs > Sgen")
# plot std error 
polygon(x=c(agg_abund, rev(agg_abund)), y= c(preds_up, rev(preds_low)), col=adjustcolor("dodgerblue", alpha=0.5), border=NA)
# Plot predictions from constrained B_0 model
lines(y=preds_c, x=agg_abund, col="gray")
polygon(x=c(agg_abund, rev(agg_abund)), y= c(preds_c_up, rev(preds_c_low)), col=adjustcolor("gray", alpha=0.5), border=NA)
# plot benchmarks
#abline(v=mod_out$Aggregate_LRPs_3phase$Estimate[mod_out$Aggregate_LRPs_3phase$param=="Agg_BM"], col="dodgerblue", lty=2)
#abline(v=mres$Estimate[mres$param=="Agg_BM"]*scale, col="gray", lty=2)

# # Get binomial regression parameters from above manual 3 phase with controls, etc.
# # plot observed data
points(y = obj$report()$N_Above_LRP/7, x= obj$report()$Agg_Abund*scale) # from manual 3 phase model fitting above

legend(x=1200000, y=0.4, col=c("dodgerblue", "gray"), lty=1, legend=c("Normal", "B_0 constrained < -2.3"))

# Plot model based on B_0 and B_1 estimates
# B_0 <- mres$Estimate[mres$param=="B_0"]
# B_1 <- mres$Estimate[mres$param=="B_1"]
# plot binomial model estimate
# curve( inv_logit(B_0 + B_1*x/scale), col="red", lty=3 , lwd=2, add=TRUE)

# Get binomial regression parameters from other model outputs and plot
# B_0 <- mod_out$Aggregate_LRPs$Estimate[mod_out$Aggregate_LRPs$param=="B_0"]
# B_1 <- mod_out$Aggregate_LRPs$Estimate[mod_out$Aggregate_LRPs$param=="B_1"]
# curve( inv_logit(B_0 + B_1*x/scale), col="orange", lty=3 , lwd=2, add=TRUE) # the one phase model results are different than the others. Does not fit data well
# B_0 <- mod_out$Aggregate_LRPs_3phase$Estimate[mod_out$Aggregate_LRPs_3phase$param=="B_0"]
# B_1 <- mod_out$Aggregate_LRPs_3phase$Estimate[mod_out$Aggregate_LRPs_3phase$param=="B_1"]
# curve( inv_logit(B_0 + B_1*x/scale), col="brown", lty=3 , lwd=2, add=TRUE) # 3 phase model gives same results 

dev.off()

# --------------------------------#
# Compare logistic fits with logit and cauchy link functions
# and with/without outlier, with generic glm function 
# --------------------------------#

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

# Plot results of TMB stock-recruit models, optional: compare ricker curves 
plot_ricker_mods(mod1 = "Aggregate_LRPs_3phase")
# Compare Ricker fits between two models
# plot_compare_mods(mod1 = models[3], mod2= "Aggregate_LRPs_3phase")

# Plot observed recruits/spawner and estimated alpha
plot_alpha_RS(mod="Aggregate_LRPs_3phase")

# Plot R/S over time
png("figures/fig_RS_by_CU.png", width=12, height=8, units="in", res=300)
ggplot(dat, aes(y=recruits/spawners, x=year)) +
  geom_point() +
  geom_path() +
  #scale_y_log10() + 
  coord_cartesian(ylim=c(0,10)) +
  geom_hline(aes(yintercept=1)) +
  facet_wrap(~CU)
dev.off()

png("figures/fig_logRS_by_CU.png", width=12, height=8, units="in", res=300)
ggplot(dat, aes(y=recruits/spawners, x=year)) +
  geom_point() +
  geom_path() +
  scale_y_log10() + 
  geom_hline(aes(yintercept=1)) +
  facet_wrap(~CU)
dev.off()

# Plot R/S over time
CUs <- unique(dat$CU)
datg <- dat
datg$group <-  ifelse(datg$CU %in% CUs[c(6,7)], "Strait of Georgia", 
                     ifelse(datg$CU %in% CUs[c(1,2,4)], "Johnstone Strait",
                            "Bute + Knight"))

# Plot R/S over time, by three major geographic groups.
# Appears that R/S track together for Strait of Georgia, Johnstone Strait. Some agreement
# between Johnstone strait and Bute/KNight, but also some very different years for those 
# two inlet CUs.
ggplot(datg, aes(y=recruits/spawners, x=year, colour= substr(CU, start=5, stop=40))) +
  geom_point() +
  geom_path(size=1.5) +
  scale_y_log10() + 
  geom_hline(aes(yintercept=1)) +
  guides(colour="none") +
  facet_wrap(~group, ncol=1)

# --------------#
# Look at linear form of ricker equation for CHum
plot_linear_ricker(mod="Aggregate_LRPs_3phase")

# Plot ricker residuals


# Plot spread in spawner abundances, chum vs. coho
# read in coho data
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
all_dat$RS <- all_dat$recruits / all_dat$spawner
# plot comparison of distributions
png(filename="figures/fig_compare_chum_coho_spawner_dist.png", width=8, height=6,units="in", res=300)
ggplot(all_dat, aes(x=spawners, colour=CU, fill=CU)) +
  geom_point(aes(y=0, x=spawners), shape=108, colour="black", size=2) +
  geom_density(alpha=0.5) + 
  scale_x_log10( breaks= c(10^(1:10))) +
  facet_wrap(~sp, ncol=1) +
  xlab(bquote("log"[10]~"(spawners)")) +
  ylab("Density") +
  theme_classic()
dev.off()

# Compare correlation among CUs
# convert long to wide data for correlations
coho_dat_w <- coho_dat %>% select(CU, spawners, year) %>% pivot_wider(names_from=CU, values_from=spawners)
chum_dat_w <- dat %>% select(CU, spawners, year) %>% pivot_wider(names_from=CU, values_from=spawners)

# Plot correlation of spawner abundances
png(filename="figures/fig_cor_spawners_coho.png", width=6, height=6, units="in", res=300)
corrplot(cor(coho_dat_w[,-1]), tl.col = "black", bg = "White",
         addCoef.col = "black", type = "upper", tl.pos="td", tl.srt=45, diag=FALSE)
#PerformanceAnalytics::chart.Correlation(coho_dat_w[,-1])
dev.off()

png(filename="figures/fig_cor_spawners_chum.png", width=10, height=10, units="in", res=300)
corrplot(cor(chum_dat_w[,-1]), tl.col = "black", bg = "White",
         addCoef.col = "black", type = "upper", tl.pos="td", tl.srt=45, diag=FALSE)
#PerformanceAnalytics::chart.Correlation(chum_dat_w[,-1])
dev.off()

# Plot correlations between recruits/spawner
# get correlations of recruits/spawner
coho_dat$RS <- coho_dat$recruits/coho_dat$spawners
coho_dat_RS_w <- coho_dat %>% select(CU, RS, year) %>% pivot_wider(names_from=CU, values_from=RS )
dat$RS <- dat$recruits/dat$spawners 
chum_dat_RS_w <- dat %>% select(CU, RS, year) %>% pivot_wider(names_from=CU, values_from=RS )
# Plot correlation of recruits/spawner
PerformanceAnalytics::chart.Correlation(coho_dat_RS_w[,-1])
PerformanceAnalytics::chart.Correlation(chum_dat_RS_w[,-1])


# -------------------------------------------#
# Save model output (to compare with Holt et al. 2018)
# -------------------------------------------#

# make a data frame with model estimates and CU names (works with any of the mod_out models)
resdf <- mod_out$Aggregate_LRPs_2phase # bind estimates with CU names
resdfw <- pivot_wider(resdf[!is.na(resdf$CU_name) , -grep("Std..Error", names(resdf))], names_from = param, values_from=Estimate) # long to wide format
resdfw$SMSY_80 <- 0.8 * resdfw$SMSY # get 80% of SMSY
#write.csv(resdfw, "output/ricker_est_to_compare_with_beta.csv", row.names=FALSE) # write to csv


res_sum <- resdfw %>% select(CU_name, A, Sgen, SMSY_80) %>% pivot_longer(cols=c(A, Sgen, SMSY_80)) # make summary table with same layout as report 
write.csv(res_sum, "output/ricker_est_to_compare.csv") # write to csv
