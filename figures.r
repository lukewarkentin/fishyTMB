# -----------------------------
# Functions to plot and save figures
#
#
#     plot_SR_by_CU - 
#     plot_compare_mods
# -----------------------------

# Plot recruits ~ spawners for each CU
plot_SR_by_CU <- function(x) {
  fig1 <- ggplot(x, aes(y=spawners, x= recruits)) + 
                   geom_point() + 
                   facet_wrap(~CU, scales="free") + 
                   theme_classic()
  ggsave("figures/fig_SR_by_CU.png", fig1)
}

#Compare ricker fits between two TMB models 
plot_compare_mods <- function(mod1, mod2) {
  CUs <- unique(dat$CU)
  png(paste0("figures/fig_compare_", mod1, "_", mod2, ".png"), height=8, width=11, units="in",res=300, pointsize=14 )
  par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
  for (i in 1:length(CUs)) {
    mod1alpha <- exp(mod_out[[mod1]][grep("logA", row.names(mod_out[[mod1]]))][i])
    mod2alpha <- exp(mod_out[[mod2]][grep("logA", row.names(mod_out[[mod2]]))][i])
    mod1beta <- exp(mod_out[[mod1]][grep("logB", row.names(mod_out[[mod1]]))][i])
    mod2beta <- exp(mod_out[[mod2]][grep("logB", row.names(mod_out[[mod2]]))][i])
    plot(dat$recruits[dat$CU==CUs[i]] ~ dat$spawners[dat$CU==CUs[i]], xlab="Spawners", ylab="Recruits")
    # plot Ricker curves:
    # R = alpha * S * exp(- beta * S)
    curve(mod1alpha * x * exp(- mod1beta * x), 
          col="black", add=TRUE)
    curve(mod2alpha * x * exp(- mod2beta * x), 
          col="dodgerblue", add=TRUE)
    # plot SMSY and Sgen
    #abline(v=res$value[grep("SMSY", names(res$value))][i], col=i)
    #abline(v=res$value[grep("Sgen", names(res$value))][i], col=i, lty=2)
    title(main = CUs[i])
  }
  plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')
  legend(x=1, y=5, legend=c(mod1, mod2), col=c("black", "dodgerblue"), lty=1)
  dev.off()
}