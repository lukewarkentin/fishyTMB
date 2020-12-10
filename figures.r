# -----------------------------
# Functions to plot and save figures
#
#
#     plot_compare_mods - plot ricker curves for two models, one panel for each stock
#     plot_ricker_mods - plot ricker curves for one model, one panel for each stock
#     plot_logistic - plot logistic model fits
# -----------------------------

# Compare ricker fits between two TMB models 
plot_compare_mods <- function(mod1, mod2) {
  CUs <- unique(dat$CU)
  png(paste0("figures/fig_compare_", mod1, "_", mod2, ".png"), height=8, width=11, units="in",res=300, pointsize=14 )
  par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
  res1 <- mod_out[[mod1]] # get data frame of results for model 1
  res2 <- mod_out[[mod2]] # get data frame of results for model 2
  for (i in 1:length(CUs)) {
    # get model 1 parameters
    mod1alpha <-     res1$Estimate[res1$param=="A"][i]
    mod1beta <-      res1$Estimate[res1$param=="B"][i] 
    # get model 2 parameters
    mod2alpha <-     res2$Estimate[res2$param=="A"][i]
    mod2beta <-      res2$Estimate[res2$param=="B"][i]
    plot(dat$recruits[dat$CU==CUs[i]] ~ dat$spawners[dat$CU==CUs[i]], xlab="Spawners", ylab="Recruits")
    # plot Ricker curves:
    # R = alpha * S * exp(- beta * S)
    curve(mod1alpha * x * exp(- mod1beta * x), 
          col="black", add=TRUE)
    curve(mod2alpha * x * exp(- mod2beta * x), 
          col="dodgerblue", add=TRUE, lty=2)
    # plot carrying capacity 1/beta
    abline(v=1/mod1beta, col="black")
    abline(v=1/mod2beta, col="dodgerblue", lty=2)
    title(main = CUs[i])
  }
  plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')
  legend(x=0, y=5, legend=c(mod1, mod2), col=c("black", "dodgerblue"), lty=c(1,2))
  dev.off()
}

# Plot ricker fits for one TMB model 
plot_ricker_mods <- function(mod1) {
  CUs <- unique(dat$CU)
  png(paste0("figures/fig_ricker_", mod1, ".png"), height=8, width=14, units="in",res=300, pointsize=14 )
  par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
  res1 <- mod_out[[mod1]] # get data frame of results for model 1
  for (i in 1:length(CUs)) {
    # get model 1 parameters
    alpha <-     res1$Estimate[res1$param=="A"][i]
    alpha_se <- res1$Std..Error[res1$param=="A"][i]
    beta <-      res1$Estimate[res1$param=="B"][i] 
    beta_se <- res1$Std..Error[res1$param=="B"][i]
    Sgen <- res1$Estimate[res1$param=="Sgen"][i]
    Sgen_se <- res1$Std..Error[res1$param=="Sgen"][i]
    SMSY <- res1$Estimate[res1$param=="SMSY"][i]
    SMSY_se <- res1$Std..Error[res1$param=="SMSY"][i]
    plot(dat$recruits[dat$CU==CUs[i]] ~ dat$spawners[dat$CU==CUs[i]], 
         xlim=c((Sgen-Sgen_se)*1.1 ,max(dat$spawners[dat$CU==CUs[i]])), xlab="Spawners", ylab="Recruits")
    # plot Ricker curves:
    # R = alpha * S * exp(- beta * S)
    # plot mean estimate curve, with upper and lower with alpha +- SE
    curve(alpha * x * exp(-beta * x), add=TRUE)
    curve((alpha + alpha_se) * x * exp(-beta * x), lty=3, add=TRUE)
    curve((alpha - alpha_se) * x * exp(-beta * x), lty=3, add=TRUE)
    # add upper and lowe with alpha +- SE and beta +- SE
    curve((alpha + alpha_se) * x * exp(-(beta-beta_se) * x), col="orange", lty=3, add=TRUE)
    curve((alpha - alpha_se) * x * exp(-(beta+beta_se) * x), col="orange", lty=3, add=TRUE)
    # Sgen with SE
    abline(v=Sgen, col="dodgerblue", lty=1)
    abline(v=Sgen+Sgen_se, col="dodgerblue", lty=3)
    abline(v=Sgen-Sgen_se, col="dodgerblue", lty=3)
    text(y=max(dat$recruits[dat$CU==CUs[i]])*0.9, x=Sgen, label="Sgen", col="dodgerblue", adj=1)
    # SMSY with SE
    abline(v=SMSY, col="firebrick", lty=1)
    abline(v=SMSY+SMSY_se, col="firebrick", lty=3)
    abline(v=SMSY-SMSY_se, col="firebrick", lty=3)
    text(y=max(dat$recruits[dat$CU==CUs[i]])*0.9, x=SMSY, label="SMSY", col="firebrick", adj=0)
    # add text of alpha
    text(y=max(dat$recruits[dat$CU==CUs[i]])*0.9, x= max(dat$spawners[dat$CU==CUs[i]]) *0.7, label=paste0("alpha = ", round(alpha,2), "\nSE = ", round(alpha_se, 2), "") )
    abline(a=0, b=1, col="gray", lty=3)
    #abline(a=0, b=alpha, col="black")
    # add 50th quantile
    # abline(v=quantile(dat$spawners[dat$CU==CUs[i]],0.5), col="purple")
    # text(y=max(dat$recruits[dat$CU==CUs[i]])*0.8, x=quantile(dat$spawners[dat$CU==CUs[i]],0.5), label="50%", col="purple", adj=0)
    title(main = CUs[i])
  }
  plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')
  dev.off()
}

# Plot observed Recruits/spawner and estimated alpha for each CU
plot_alpha_RS <- function(mod) {
  CUs <- unique(dat$CU)
  png(paste0("figures/fig_alpha_RS_", mod, ".png"), height=8, width=14, units="in",res=300, pointsize=14 )
  par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
  res1 <- mod_out[[mod]] # get data frame of results for model 1
  for (i in 1:length(CUs)) {
    # get model 1 parameters
    alpha <-     res1$Estimate[res1$param=="A"][i]
    alpha_se <-  res1$Std..Error[res1$param=="A"][i]
    # plot observed recruits/spawner
    plot(density(dat$recruits[dat$CU==CUs[i]] / dat$spawners[dat$CU==CUs[i]]), type="l", main="",
          xlab="log(recruits/spawner)")
    abline(v=alpha, col="dodgerblue")
    abline(v=alpha+alpha_se, col="dodgerblue", lty=3)
    abline(v=alpha-alpha_se, col="dodgerblue", lty=3)
    abline(v=1, col="orange", lty=3)
    title(main = CUs[i])
  }
  plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')
  legend(x=1,y=5,col=c("black", "dodgerblue"), lty=c(1,1), legend=c("Observed R/S", "alpha from Ricker"))
  dev.off()
}

# Plot observed log(Recruits/spawner) and against spawners, with alpha and beta estimates
plot_linear_ricker <- function(mod) {
  CUs <- unique(dat$CU)
  png(paste0("figures/fig_linear_ricker_", mod, ".png"), height=8, width=14, units="in",res=300, pointsize=14 )
  par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
  res1 <- mod_out[[mod]] # get data frame of results for model 1
  for (i in 1:length(CUs)) {
    # get model 1 parameters
    alpha <-     res1$Estimate[res1$param=="A"][i]
    beta <- res1$Estimate[res1$param=="B"][i]
    SMSY <- res1$Estimate[res1$param=="SMSY"][i]
    Sgen <- res1$Estimate[res1$param=="Sgen"][i]
    # plot observed recruits/spawner
    plot( log(dat$recruits[dat$CU==CUs[i]] / dat$spawners[dat$CU==CUs[i]]) ~ dat$spawners[dat$CU==CUs[i]], type="p", main="",
         xlab="spawners", ylab="log(recruits/spawner)")
    abline(a=log(alpha), b=-beta, col="dodgerblue")
    abline(h=0, col="gray", lty=3)
    abline(v=Sgen, col="orange", lty=2)
    abline(v=SMSY, col="purple", lty=2)
    title(main = CUs[i])
  }
  plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')
  legend(x=1,y=5, lty=2, col=c("orange", "purple"), legend=c("Sgen", "SMSY"))
  dev.off()
}

