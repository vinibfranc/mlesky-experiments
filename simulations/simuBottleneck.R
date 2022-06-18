#rm(list=ls())

# load(file="simuBottle.RData")
source("simulations/commonFunctions.R")
start <- Sys.time()

#bottleFun=function(x){if (x<2005 || x>2010) 10 else 1}
bottleFun=function(x){
  res = vector()
  for(i in 1:length(x)) {
    #print(x[[i]])
    if(x[i]<2005 || x[i]>2010)
      res[i]=10
    else
      res[i]=1
  }
  return(res)
}

t200_bottle_trees = simulate_ntrees(sampleDates_200, bottleFun)
t50_bottle_trees = simulate_ntrees(sampleDates_50, bottleFun)
t20_bottle_trees = simulate_ntrees(sampleDates_20, bottleFun)

t200_bottle_m1 = get_nsim_estimates(t200_bottle_trees, sampleDates_200, bottleFun, 1, "bottle/bottle_t200_m1.png", "bottle/bottle_t200_m1_mae.png", "bottle/bottle_t200_m1_rmse.png")
t200_bottle_m2 = get_nsim_estimates(t200_bottle_trees, sampleDates_200, bottleFun, 2, "bottle/bottle_t200_m2.png", "bottle/bottle_t200_m2_mae.png", "bottle/bottle_t200_m2_rmse.png")
t200_bottle_m3 = get_nsim_estimates(t200_bottle_trees, sampleDates_200, bottleFun, 3, "bottle/bottle_t200_m3.png", "bottle/bottle_t200_m3_mae.png", "bottle/bottle_t200_m3_rmse.png")

t50_bottle_m1 = get_nsim_estimates(t50_bottle_trees, sampleDates_50, bottleFun, 1, "bottle/bottle_t50_m1.png", "bottle/bottle_t50_m1_mae.png", "bottle/bottle_t50_m1_rmse.png")
t50_bottle_m2 = get_nsim_estimates(t50_bottle_trees, sampleDates_50, bottleFun, 2, "bottle/bottle_t50_m2.png", "bottle/bottle_t50_m2_mae.png", "bottle/bottle_t50_m2_rmse.png")
t50_bottle_m3 = get_nsim_estimates(t50_bottle_trees, sampleDates_50, bottleFun, 3, "bottle/bottle_t50_m3.png", "bottle/bottle_t50_m3_mae.png", "bottle/bottle_t50_m3_rmse.png")

t20_bottle_m1 = get_nsim_estimates(t20_bottle_trees, sampleDates_20, bottleFun, 1, "bottle/bottle_t20_m1.png", "bottle/bottle_t20_m1_mae.png", "bottle/bottle_t20_m1_rmse.png")
t20_bottle_m2 = get_nsim_estimates(t20_bottle_trees, sampleDates_20, bottleFun, 2, "bottle/bottle_t20_m2.png", "bottle/bottle_t20_m2_mae.png", "bottle/bottle_t20_m2_rmse.png")
t20_bottle_m3 = get_nsim_estimates(t20_bottle_trees, sampleDates_20, bottleFun, 3, "bottle/bottle_t20_m3.png", "bottle/bottle_t20_m3_mae.png", "bottle/bottle_t20_m3_rmse.png")

# compare different models for the same function and sample size
t200_bottle_models = compare_cov_prob_models(t200_bottle_m1[[7]], t200_bottle_m2[[7]], t200_bottle_m3[[7]], t200_bottle_m1[[8]], "comp_models/cov_prob/bottle_t200.png")
t50_bottle_models = compare_cov_prob_models(t50_bottle_m1[[7]], t50_bottle_m2[[7]], t50_bottle_m3[[7]], t50_bottle_m1[[8]], "comp_models/cov_prob/bottle_t50.png")
t20_bottle_models = compare_cov_prob_models(t20_bottle_m1[[7]], t20_bottle_m2[[7]], t20_bottle_m3[[7]], t20_bottle_m1[[8]], "comp_models/cov_prob/bottle_t20.png")

ggarrange(t200_bottle_models, t50_bottle_models, t20_bottle_models, nrow=3, ncol=1, labels=c("n=200","n=50","n=20"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("coverage_plots/comp_models/bottle_models.pdf", units="in", width=9, height=12, dpi=600)

# compare different sample sizes for same function and model
tx_bottle_skysigma = compare_same_model_diff_samp_size(t200_bottle_m1[[7]], t50_bottle_m1[[7]], t20_bottle_m1[[7]], t200_bottle_m1[[8]], "comp_models/cov_prob/bottle_diff_samp_skysigma.png")
tx_bottle_skygrid = compare_same_model_diff_samp_size(t200_bottle_m2[[7]], t50_bottle_m2[[7]], t20_bottle_m2[[7]], t200_bottle_m2[[8]], "comp_models/cov_prob/bottle_diff_samp_skygrid.png")
tx_bottle_skygrowth = compare_same_model_diff_samp_size(t200_bottle_m3[[7]], t50_bottle_m3[[7]], t20_bottle_m3[[7]], t200_bottle_m3[[8]], "comp_models/cov_prob/bottle_diff_samp_skygrowth.png")

ggarrange(tx_bottle_skysigma, tx_bottle_skygrid, tx_bottle_skygrowth, nrow=3, ncol=1, labels=c("Skysigma","Skygrid","Skygrowth"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("coverage_plots/comp_models/bottle_diff_samp_sizes.pdf", units="in", width=9, height=12, dpi=600)

# compare error considering different models for the same function and sample size
t200_bottle_models_err = compare_err_models(t200_bottle_m1[[9]], t200_bottle_m2[[9]], t200_bottle_m3[[9]], t200_bottle_m1[[10]], t200_bottle_m2[[10]], t200_bottle_m3[[10]], t200_bottle_m1[[8]], "comp_models/mae/bottle_t200.png", "comp_models/rmse/bottle_t200.png")
t50_bottle_models_err = compare_err_models(t50_bottle_m1[[9]], t50_bottle_m2[[9]], t50_bottle_m3[[9]], t50_bottle_m1[[10]], t50_bottle_m2[[10]], t50_bottle_m3[[10]], t50_bottle_m1[[8]], "comp_models/mae/bottle_t50.png", "comp_models/rmse/bottle_t50.png")
t20_bottle_models_err = compare_err_models(t20_bottle_m1[[9]], t20_bottle_m2[[9]], t20_bottle_m3[[9]], t20_bottle_m1[[10]], t20_bottle_m2[[10]], t20_bottle_m3[[10]], t20_bottle_m1[[8]], "comp_models/mae/bottle_t20.png", "comp_models/rmse/bottle_t20.png")

ggarrange(t200_bottle_models_err[[1]], t50_bottle_models_err[[1]], t20_bottle_models_err[[1]], nrow=3, ncol=1, labels=c("n=200","n=50","n=20"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("error_plots/comp_models/mae/bottle_models.pdf", units="in", width=9, height=12, dpi=600)

ggarrange(t200_bottle_models_err[[2]], t50_bottle_models_err[[2]], t20_bottle_models_err[[2]], nrow=3, ncol=1, labels=c("n=200","n=50","n=20"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("error_plots/comp_models/rmse/bottle_models.pdf", units="in", width=9, height=12, dpi=600)

# compare error considering different sample sizes for same function and model
tx_bottle_skysigma_err = compare_err_same_model_diff_samp_size(t200_bottle_m1[[9]], t50_bottle_m1[[9]], t20_bottle_m1[[9]], t200_bottle_m1[[10]], t50_bottle_m1[[10]], t20_bottle_m1[[10]], t200_bottle_m1[[8]], "comp_models/mae/bottle_diff_samp_skysigma.png", "comp_models/rmse/bottle_diff_samp_skysigma.png")
tx_bottle_skygrid_err = compare_err_same_model_diff_samp_size(t200_bottle_m2[[9]], t50_bottle_m2[[9]], t20_bottle_m2[[9]], t200_bottle_m2[[10]], t50_bottle_m2[[10]], t20_bottle_m2[[10]], t200_bottle_m2[[8]], "comp_models/mae/bottle_diff_samp_skygrid.png", "comp_models/rmse/bottle_diff_samp_skygrid.png")
tx_bottle_skygrowth_err = compare_err_same_model_diff_samp_size(t200_bottle_m3[[9]], t50_bottle_m3[[9]], t20_bottle_m3[[9]], t200_bottle_m3[[10]], t50_bottle_m3[[10]], t20_bottle_m3[[10]], t200_bottle_m3[[8]], "comp_models/mae/bottle_diff_samp_skygrowth.png", "comp_models/rmse/bottle_diff_samp_skygrowth.png")

ggarrange(tx_bottle_skysigma_err[[1]], tx_bottle_skygrid_err[[1]], tx_bottle_skygrowth_err[[1]], nrow=3, ncol=1, labels=c("Skysigma","Skygrid","Skygrowth"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("error_plots/comp_models/mae/bottle_diff_samp_sizes.pdf", units="in", width=9, height=12, dpi=600)

ggarrange(tx_bottle_skysigma_err[[2]], tx_bottle_skygrid_err[[2]], tx_bottle_skygrowth_err[[2]], nrow=3, ncol=1, labels=c("Skysigma","Skygrid","Skygrowth"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("error_plots/comp_models/rmse/bottle_diff_samp_sizes.pdf", units="in", width=9, height=12, dpi=600)

end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
print(paste("Total time elapsed: ",total_time,"mins"))

save.image(file="simuBottle.RData")
