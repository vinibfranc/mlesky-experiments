#rm(list=ls())

# load(file="simuSinus.RData")
source("simulations/commonFunctions.R")
start <- Sys.time()

sinusFun=function(x){sin(x)*10+12}

t200_sinus_trees = simulate_ntrees(sampleDates_200, sinusFun)
t50_sinus_trees = simulate_ntrees(sampleDates_50, sinusFun)
t20_sinus_trees = simulate_ntrees(sampleDates_20, sinusFun)

t200_sinus_m1 = get_nsim_estimates(t200_sinus_trees, sampleDates_200, sinusFun, 1, "sinus/sinus_t200_m1.png")
t200_sinus_m2 = get_nsim_estimates(t200_sinus_trees, sampleDates_200, sinusFun, 2, "sinus/sinus_t200_m2.png")
t200_sinus_m3 = get_nsim_estimates(t200_sinus_trees, sampleDates_200, sinusFun, 3, "sinus/sinus_t200_m3.png")

t50_sinus_m1 = get_nsim_estimates(t50_sinus_trees, sampleDates_50, sinusFun, 1, "sinus/sinus_t50_m1.png")
t50_sinus_m2 = get_nsim_estimates(t50_sinus_trees, sampleDates_50, sinusFun, 2, "sinus/sinus_t50_m2.png")
t50_sinus_m3 = get_nsim_estimates(t50_sinus_trees, sampleDates_50, sinusFun, 3, "sinus/sinus_t50_m3.png")

t20_sinus_m1 = get_nsim_estimates(t20_sinus_trees, sampleDates_20, sinusFun, 1, "sinus/sinus_t20_m1.png")
t20_sinus_m2 = get_nsim_estimates(t20_sinus_trees, sampleDates_20, sinusFun, 2, "sinus/sinus_t20_m2.png")
t20_sinus_m3 = get_nsim_estimates(t20_sinus_trees, sampleDates_20, sinusFun, 3, "sinus/sinus_t20_m3.png")

# compare different models for the same function and sample size
t200_sinus_models = compare_cov_prob_models(t200_sinus_m1[[7]], t200_sinus_m2[[7]], t200_sinus_m3[[7]], t200_sinus_m1[[8]], "comp_models/sinus_t200.png")
t50_sinus_models = compare_cov_prob_models(t50_sinus_m1[[7]], t50_sinus_m2[[7]], t50_sinus_m3[[7]], t50_sinus_m1[[8]], "comp_models/sinus_t50.png")
t20_sinus_models = compare_cov_prob_models(t20_sinus_m1[[7]], t20_sinus_m2[[7]], t20_sinus_m3[[7]], t20_sinus_m1[[8]], "comp_models/sinus_t20.png")

ggarrange(t200_sinus_models, t50_sinus_models, t20_sinus_models, nrow=3, ncol=1, labels=c("n=200","n=50","n=20"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("coverage_plots/comp_models/sinus_models.pdf", units="in", width=9, height=12, dpi=600)

# compare different sample sizes for same function and model
tx_sinus_skysigma = compare_same_model_diff_samp_size(t200_sinus_m1[[7]], t50_sinus_m1[[7]], t20_sinus_m1[[7]], t200_sinus_m1[[8]], "comp_models/sinus_diff_samp_skysigma.png")
tx_sinus_skygrid = compare_same_model_diff_samp_size(t200_sinus_m2[[7]], t50_sinus_m2[[7]], t20_sinus_m2[[7]], t200_sinus_m2[[8]], "comp_models/sinus_diff_samp_skygrid.png")
tx_sinus_skygrowth = compare_same_model_diff_samp_size(t200_sinus_m3[[7]], t50_sinus_m3[[7]], t20_sinus_m3[[7]], t200_sinus_m3[[8]], "comp_models/sinus_diff_samp_skygrowth.png")

ggarrange(tx_sinus_skysigma, tx_sinus_skygrid, tx_sinus_skygrowth, nrow=3, ncol=1, labels=c("Skysigma","Skygrid","Skygrowth"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("coverage_plots/comp_models/sinus_diff_samp_sizes.pdf", units="in", width=9, height=12, dpi=600)

end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
print(paste("Total time elapsed: ",total_time,"mins"))

save.image(file="simuSinus.RData")