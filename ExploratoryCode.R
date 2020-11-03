#library(R2jags)
library(TMB)
library(tidyverse)
source("Code/Functions.R")

# Troubleshooting problem data sets that TMB model gives weird LRP's for


Sim_Out <- readRDS( "Simulations/DataOut/ProblemData.rds")
Sim_Params <- Sim_Out$Sim_Params
Sim_Data <- Sim_Out$Sim_Data

LRP_Summ <- left_join(Sim_Data, Sim_Params) %>% group_by(Year) %>% summarise(Total_ETS = sum(ets),
                                      CUs_Above_LRP = sum(ets > True_SGen),
                                      n_CUs = n(),
                                      Prop_Above_LRP = sum(ets > True_SGen)/n() )
Prop_Mod = glm(cbind(CUs_Above_LRP, n_CUs - CUs_Above_LRP) ~ Total_ETS, family = binomial, data = LRP_Summ)
GLM_LRP <- (log(0.95/(1-0.95)) - Prop_Mod$coefficients[1])/ Prop_Mod$coefficients[2]

# xx values to predict from
xx <- data.frame(Total_ETS = seq(0, max(LRP_Summ$Total_ETS*1.25), by=(max(LRP_Summ$Total_ETS*1.25)/1000)))
Prop_yy <- predict(Prop_Mod, xx, type="response")
plot(xx$Total_ETS, Prop_yy, col="red", ylim=c(0,1), type="o", pch=19, xlab="", ylab="", cex.axis=1.2)
mtext(side=1, "Aggregate Escapement", line=3,las=0)
mtext(side = 2, "Proportion CUs > LRP", line=3,las=0)
points(LRP_Summ$Total_ETS, LRP_Summ$Prop_Above_LRP, pch=19)   
abline(v=GLM_LRP, col = "orange", lwd=2, lty=2)

# compare this to TMB estimate
Data <- Sim_Data
data <- list()
data$S <- Data$ets/Scale #Data$Scaled_ETS
#data$S <- Data$Scaled_ETS
data$logR <- log(Data$rec/Scale)
#data$logR <- log(Data$Scaled_R)
data$stk <- as.numeric(Data$Stock_Name)-1
N_Stocks <- length(unique(Data$Stock_Name))
data$N_Stks <- N_Stocks
data$yr <- Data$Year
Mod_Yrs <- Data$Year
data$Mod_Yr_0 <- min(Mod_Yrs)
data$Mod_Yr_n <- max(Mod_Yrs)
#data$Scale <- Scales

param <- list()
param$logA <- rep(0.5, N_Stocks)
param$logB <- log(1/((Data %>% group_by(Stock_Name) %>%  summarise(x=quantile(ets, 0.75)))$x/Scale)) 
#param$logB <- log(1/((Data %>% group_by(Stock_Name) %>%  summarise(x=quantile(Scaled_ETS, 0.75)))$x))
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <-  log((Data %>% group_by(Stock_Name) %>%  summarise(x=quantile(ets, 0.25)))$x/Scale) 
#param$logSgen <-  log((Data %>% group_by(Stock_Name) %>%  summarise(x=quantile(Scaled_ETS, 0.25)))$x)
param$B_0 <- 2
param$B_1 <- 0.1

if(Mod == "Aggregate_LRPs_Hier"){
  param$logMuA <- 3
  param$logSigmaA <- -2
}

# Phase 1 estimate SR params
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
if(Mod == "Aggregate_LRPs"){
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, map=map)
} else if(Mod == "Aggregate_LRPs_Hier"){
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)
}
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj))

# pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
upper = c(rep(Inf, 3*N_Stocks), log(SMSYs), Inf, Inf )
if (Mod ==  "Aggregate_LRPs_Hier"){
  upper <- c(upper, c(Inf, Inf))
}

#Phase 2 get Sgen, SMSY etc.
map2 <- list(B_0=factor(NA), B_1=factor(NA))
if(Mod == "Aggregate_LRPs"){
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, map=map2)
} else if(Mod == "Aggregate_LRPs_Hier"){
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA", map=map2)
}
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper )
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
summary(sdreport(obj))

# Phase 3 fit logistic model
# Hold other estimates constant
if(Mod == "Aggregate_LRPs"){
  obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE)
} else if(Mod == "Aggregate_LRPs_Hier"){
  map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
  obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE, map=map3)
}
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper)
pl3 <- obj$env$parList(opt$par)
summary(sdreport(obj))

 out<- obj$report()
 Prop_Mod2 = glm(cbind(N_Above_LRP, N_Stocks - N_Above_LRP) ~ Agg_Abund, family = binomial, data =out)
 # Get LRP (calculated as aggregate ETS number associated with 95% probability of all > LRP)?
 LRP2 <- (log(0.95/(1-0.95)) - Prop_Mod2$coefficients[1])/ Prop_Mod2$coefficients[2]
 # GLM actually gives same value -- so not TMB issue in this case
 
 points(out$Agg_Abund*10000, out$N_Above_LRP/10, col="blue") # output looks similar but not exactly the same
 xx <- data.frame(Agg_Abund = seq(0, max(out$Agg_Abund*1.25), by=(max(out$Agg_Abund*1.25)/1000)))
 Prop_yy2 <- predict(Prop_Mod2, xx, type="response")
 
 lines(xx$Agg_Abund*10000, Prop_yy2, col="green")

All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)

All_Ests1 <- data.frame(summary(sdreport(obj)))
All_Ests1$Param <- row.names(All_Ests1)
yy1 <- inv_logit( All_Ests1$Estimate[All_Ests1$Param == "B_0"] + 
                    All_Ests1$Estimate[All_Ests1$Param == "B_1"]*xx)


Fits <- list()
Fits$Stock_Ests <- data.frame(Stock_Name = as.factor(1:N_Stocks), Est_A = exp(pl3$logA), Est_B = exp(pl3$logB)/Scale, Sig_Est = exp(pl3$logSigma),
                              SMSY_Est = All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]*Scale, 
                              Sgen_Est = exp(pl3$logSgen)*Scale)
Fits$Agg_Ests <- All_Ests[All_Ests$Param == "Agg_BM",]
Fits$Agg_Ests <-  Fits$Agg_Ests[,c("Estimate", "Std..Error")]*Scale





















compile("Code/Basic_Ricker.cpp")
dyn.load(dynlib("Code/Basic_Ricker"))
compile("Code/Hier_Ricker.cpp")
dyn.load(dynlib("Code/Hier_Ricker"))

# already simulated set of 7 stocks with alphas coming from normal dist (mean=4, sd=1)

Sim_Data <- readRDS("Code/Sim_Data.RDS")

# simulate stocks with higher uncert
N_Stocks = 7
true_a <- abs(rnorm(N_Stocks, 4, 2 ))
true_b <- 1/runif(N_Stocks, 2000, 15000)
leng <- 75

Sim_Data <- data.frame(Year = numeric(), Stock_Name = numeric(), Escape = numeric(), Rec = numeric(), 
                       True_A = numeric(), True_B = numeric(), True_SGen = numeric())

for(i in 1:N_Stocks){
  Out <- fake_SR_data(leng=leng, true_a = true_a[i], true_b = true_b[i], sigma = 0.6,
                      hr_mean = 0.5, hr_sd = 0.2, hr_max = 0.85, init_min = 1000, init_max = 3000)
  New_Rows <- data.frame(Year = 1:leng, Stock_Name = rep(i, leng), ets = round(Out$esc), 
                         rec = round(Out$rec), True_A = rep(Out$true_a, leng), 
                         True_B = rep(Out$true_b, leng), True_SGen = rep(Out$SGen, leng))
  Sim_Data <- rbind(Sim_Data, New_Rows)
}

ggplot(Sim_Data, aes(x=ets, y=rec)) +geom_point()+ facet_wrap(~Stock_Name)

Comp_Hier(Sim_Data, fname="SigA2_SigR0_6")