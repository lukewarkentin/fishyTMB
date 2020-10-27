# -----------------------------
# Code to plot and save figures
# -----------------------------


# Plot recruits ~ spawners for each CU
plot_SR_by_CU <- function(x) {
  fig1 <- ggplot(x, aes(y=spawners, x= recruits)) + 
                   geom_point() + 
                   facet_wrap(~CU, scales="free") + 
                   theme_classic()
  ggsave("figures/fig_SR_by_CU.png", fig1)
}

