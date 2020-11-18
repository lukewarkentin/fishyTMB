# -----------------------------
# Functions to plot and save figures
#
#
#     plot_compare_mods - plot ricker curves for two models, one panel for each stock
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