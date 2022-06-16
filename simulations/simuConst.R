#rm(list=ls())

# load(file="simuConst.RData")
source("simulations/commonFunctions.R")
start <- Sys.time()

constFun=function(x){20}

t200_const_trees = simulate_ntrees(sampleDates_200, constFun)
t50_const_trees = simulate_ntrees(sampleDates_50, constFun)
t20_const_trees = simulate_ntrees(sampleDates_20, constFun)

t200_const_m1 = get_nsim_estimates(t200_const_trees, sampleDates_200, constFun, 1, "const/const_t200_m1.png")
t200_const_m2 = get_nsim_estimates(t200_const_trees, sampleDates_200, constFun, 2, "const/const_t200_m2.png")
t200_const_m3 = get_nsim_estimates(t200_const_trees, sampleDates_200, constFun, 3, "const/const_t200_m3.png")

t50_const_m1 = get_nsim_estimates(t50_const_trees, sampleDates_50, constFun, 1, "const/const_t50_m1.png")
t50_const_m2 = get_nsim_estimates(t50_const_trees, sampleDates_50, constFun, 2, "const/const_t50_m2.png")
t50_const_m3 = get_nsim_estimates(t50_const_trees, sampleDates_50, constFun, 3, "const/const_t50_m3.png")

t20_const_m1 = get_nsim_estimates(t20_const_trees, sampleDates_20, constFun, 1, "const/const_t20_m1.png")
t20_const_m2 = get_nsim_estimates(t20_const_trees, sampleDates_20, constFun, 2, "const/const_t20_m2.png")
t20_const_m3 = get_nsim_estimates(t20_const_trees, sampleDates_20, constFun, 3, "const/const_t20_m3.png")

# compare different models for the same function and sample size
t200_const_models = compare_cov_prob_models(t200_const_m1[[7]], t200_const_m2[[7]], t200_const_m3[[7]], t200_const_m1[[8]], "comp_models/const_t200.png")
t50_const_models = compare_cov_prob_models(t50_const_m1[[7]], t50_const_m2[[7]], t50_const_m3[[7]], t50_const_m1[[8]], "comp_models/const_t50.png")
t20_const_models = compare_cov_prob_models(t20_const_m1[[7]], t20_const_m2[[7]], t20_const_m3[[7]], t20_const_m1[[8]], "comp_models/const_t20.png")

ggarrange(t200_const_models, t50_const_models, t20_const_models, nrow=3, ncol=1, labels=c("n=200","n=50","n=20"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("coverage_plots/comp_models/const_models.pdf", units="in", width=9, height=12, dpi=600)

# compare different sample sizes for same function and model
tx_const_skysigma = compare_same_model_diff_samp_size(t200_const_m1[[7]], t50_const_m1[[7]], t20_const_m1[[7]], t200_const_m1[[8]], "comp_models/const_diff_samp_skysigma.png")
tx_const_skygrid = compare_same_model_diff_samp_size(t200_const_m2[[7]], t50_const_m2[[7]], t20_const_m2[[7]], t200_const_m2[[8]], "comp_models/const_diff_samp_skygrid.png")
tx_const_skygrowth = compare_same_model_diff_samp_size(t200_const_m3[[7]], t50_const_m3[[7]], t20_const_m3[[7]], t200_const_m3[[8]], "comp_models/const_diff_samp_skygrowth.png")

ggarrange(tx_const_skysigma, tx_const_skygrid, tx_const_skygrowth, nrow=3, ncol=1, labels=c("Skysigma","Skygrid","Skygrowth"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("coverage_plots/comp_models/const_diff_samp_sizes.pdf", units="in", width=9, height=12, dpi=600)

end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
print(paste("Total time elapsed: ",total_time,"mins"))

save.image(file="simuConst.RData")