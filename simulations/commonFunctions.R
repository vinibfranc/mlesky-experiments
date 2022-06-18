library(devtools)
# Install mlesky updated version including parboot modifications v0.1.4 (27 May 2022) 
#devtools::install_github("emvolz-phylodynamics/mlesky")
library(mlesky)
library(ape)
#install.packages("pbmcapply")
library(pbmcapply)
library(dplyr)
library(ggplot2)
library(ggpubr)
#library(data.table)
library(reshape2)
#install.packages("Metrics")
#library(Metrics)
library(stringr)

ncpu = 7
nsim = 500
scaleFUN = function(x) sprintf("%.2f", x) # format to 2 decimal places

set.seed(3)
#set.seed(425672)
sampleDates_200 =seq(2000,2020,0.1)
sampleDates_50 =seq(2000,2020,0.4)
sampleDates_20 =seq(2004,2014,0.5)
# sampleDates_20 =seq(2005,2010,0.25)

# Simulate coalescent trees using the simCoal method implemented in mlesky
simulate_ntrees <- function(sampDates, alphaFun) {
  tres = pbmclapply(1:nsim, function(i) {
    simCoal(sampDates, alphaFun, alphaMin=0.1)
  }, mc.cores = ncpu)
  return(tres)
}

# Get Ne estimates on common time axis and plot coverage probabilities over time
# Return a list contaning 8 items:
# (i) pboot_ne_ind: Ne estimates for each tree without getting common time axis
# (ii) pboot_ne_adj_df: Ne estimates using common time axis
# (iii) cov_prob: all Ne estimates and true value and its coverages on dataframe
# (iv) cov_prob_matrix: coverage probability matrix (rows = time axis elements; columns = simulation indexes)
# (v) cov_prob_avg: coverage probability over entire time axis
# (vi) cov_prob_over_time: coverage probability for each common time axis point over different simulations (rowMeans output)
# (vii) cov_prob_over_time_df: same as (vi) but on data.frame format to allow easy plotting using ggplot2
# (viii) common_time_ax: common time axis to use for plots
# (ix) mean_abs_error_df: Mean absolute error (MAE) for each simulation index
# (x) rmse_df: Root mean squared error (RMSE) for each simulation index
# Also plot the coverage probability over time to `coverage_plots` folder that is appended to the supplied `out_path`
# and print most of these returned values on the screen in the end of execution
get_nsim_estimates <- function(sim_trees, sampDates, alphaFun, model, out_path_cov_prob, out_path_error, out_path_rmse) {
  print("Fit")
  sts = sampDates; names(sts) = sim_trees[[1]]$tip.label
  fit = pbmclapply(X=sim_trees, FUN=mlskygrid, sampleTimes=sts, res=NULL,tau=NULL,tau_lower = 0.001,tau_upper = 1000,model = model,ncpu=ncpu)
  pboot = pboot_ne = list()
  for(i in 1:length(sim_trees)){
    print(paste("===",i,"==="))
    #print(fit[[i]]$time)
    #print(fit[[i]]$ne)
    if(length(fit[[i]]$time) <= 2) {
      next
    }
    pboot[[i]] = parboot(fit[[i]], nrep=200, ncpu=ncpu)
    pboot_ne[[i]] = pboot[[i]]["ne_ci"]
    pboot_ne[[i]]$nelb = pboot[[i]]$ne_ci[,1]
    pboot_ne[[i]]$ne = pboot[[i]]$ne_ci[,2]
    pboot_ne[[i]]$neub = pboot[[i]]$ne_ci[,3]
    pboot_ne[[i]]$time = pboot[[i]]$time; pboot_ne[[i]]$time = as.data.frame(pboot_ne[[i]]$time)
    #pboot_ne[[i]]$true_ne = alphaFun(pboot_ne[[i]]$time); #pboot_ne[[i]]$true_ne = as.data.frame(pboot_ne[[i]]$true_ne)
    pboot_ne[[i]] = cbind(pboot_ne[[i]]$time, pboot_ne[[i]]$ne_ci) #pboot_ne[[i]]$true_ne
    colnames(pboot_ne[[i]]) = c("time","nelb","est_ne","neub") #"true_ne"
  }
  # Get estimates for each tree without getting common time axis
  pboot_ne_ind = bind_rows(pboot_ne, .id = 'n_sim')
  print("Estimates without common time axis:")
  print(pboot_ne_ind)

  # Get first most recent coalescent time and last less recent coalescent time among all simulated trees to plot out Ne values on common time axis
  first_coal_sims = pboot_ne_ind %>% group_by(n_sim) %>% filter(row_number()==1)
  mr_first_coal_time = max(first_coal_sims$time)
  
  last_coal_sims = pboot_ne_ind %>% group_by(n_sim) %>% filter(row_number()==n())
  oldest_last_coal_time = min(last_coal_sims$time)
  
  # Get common time axis for all simulations
  common_time_ax = approx(pboot_ne_ind$time, pboot_ne_ind$est_ne, xout=seq(mr_first_coal_time, oldest_last_coal_time, length.out=nrow(pboot_ne_ind)/length(sim_trees)), rule=2)$x
  print("Common time axis:")
  print(common_time_ax)
  
  # Split df by simulation index again and compute approximate Ne estimates based on common time axis
  pboot_ne_ind_common = split(pboot_ne_ind, pboot_ne_ind$n_sim)
  time_adj = true_ne_adj = ne_est_adj = nelb_adj = neub_adj = pboot_ne_adj = list()
  for(i in 1:length(pboot_ne_ind_common)) { 
    time_adj[[i]] = common_time_ax; #time_adj[[i]] = as.data.frame(time_adj[[i]])
    true_ne_adj[[i]] = alphaFun(time_adj[[i]]); #true_ne_adj[[i]] = as.data.frame(true_ne_adj[[i]])
    ne_est_adj[[i]] = approx(pboot_ne_ind_common[[i]]$time, pboot_ne_ind_common[[i]]$est_ne, xout=common_time_ax, rule=2)$y; ne_est_adj[[i]] = as.data.frame(ne_est_adj[[i]])
    nelb_adj[[i]] = approx(pboot_ne_ind_common[[i]]$time, pboot_ne_ind_common[[i]]$nelb, xout=common_time_ax, rule=2)$y; nelb_adj[[i]] = as.data.frame(nelb_adj[[i]])
    neub_adj[[i]] = approx(pboot_ne_ind_common[[i]]$time, pboot_ne_ind_common[[i]]$neub, xout=common_time_ax, rule=2)$y; neub_adj[[i]] = as.data.frame(neub_adj[[i]])
    # mean_abs_error[[i]] = sum(abs(true_ne_adj[[i]] - ne_est_adj[[i]]))/length(true_ne_adj[[i]])
    # mae_metrics[[i]] = mae(actual = true_ne_adj[[i]], predicted = ne_est_adj[[i]])
    # rmse_val[[i]] = sqrt(mean((true_ne_adj[[i]] - ne_est_adj[[i]])^2))
    # rmse_metrics[[i]] = rmse(true_ne_adj[[i]], ne_est_adj[[i]])
    pboot_ne_adj[[i]] = cbind(time_adj[[i]], true_ne_adj[[i]], nelb_adj[[i]], ne_est_adj[[i]], neub_adj[[i]])
    colnames(pboot_ne_adj[[i]]) = c("time","true_ne", "nelb", "est_ne", "neub")
  }
  # Join again now with common time axis
  pboot_ne_adj_df = bind_rows(pboot_ne_adj, .id = 'n_sim')
  print("Estimates using common time axis:")
  print(pboot_ne_adj_df)
  
  # Calculate Mean Absolute Error (MEA) and Root Mean Square Error (RMSE)
  mean_abs_error_df = pboot_ne_adj_df %>% group_by(n_sim) %>% summarise(mean_abs_error = sum(abs(true_ne - est_ne))/length(common_time_ax))
  lvls_mae <- str_sort(unique(mean_abs_error_df$n_sim), numeric = TRUE)
  mean_abs_error_df$n_sim <- factor(mean_abs_error_df$n_sim, levels = lvls_mae)
  #mean_abs_error_df$n_sim = factor(mean_abs_error_df$n_sim, levels=unique(mean_abs_error_df$n_sim))
  rmse_df = pboot_ne_adj_df %>% group_by(n_sim) %>% summarise(rmse_val = sqrt(mean((true_ne - est_ne)^2)))
  lvls_rmse <- str_sort(unique(rmse_df$n_sim), numeric = TRUE)
  rmse_df$n_sim <- factor(rmse_df$n_sim, levels = lvls_rmse)
  #rmse_df$n_sim = factor(rmse_df$n_sim, levels=unique(rmse_df$n_sim))
  # Got same results from Metrics package functions
  #mae_metrics = pboot_ne_adj_df %>% group_by(n_sim) %>% summarise(mae(actual = true_ne, predicted = est_ne))
  #rmse_metrics = pboot_ne_adj_df %>% group_by(n_sim) %>% summarise(rmse(actual = true_ne, predicted = est_ne))
  # Summary error and RMSE over all simulations
  mean_abs_error_avg = mean(mean_abs_error_df$mean_abs_error)
  rmse_avg = mean(rmse_df$rmse_val)
  
  # Calculate coverage probability (whether the true Ne value is covered [=1] or not [=0] by the estimated confidence intervals)
  cov_prob = pboot_ne_adj_df %>% mutate(cover = ifelse((true_ne >= nelb & true_ne <= neub), 1, 0))
  print("Coverage probability full df:")
  print(cov_prob)
  
  cov_prob_min = cov_prob %>% select(n_sim, time, cover)
  cov_prob_min = dcast(cov_prob_min,time~n_sim, value.var = "cover")
  print(cov_prob_min)
  cov_prob_matrix = as.matrix(cov_prob_min)
  rownames(cov_prob_matrix) = cov_prob_matrix[,1]
  cov_prob_matrix = cov_prob_matrix[,-1]
  print("Coverage probability matrix:")
  print(cov_prob_matrix)
  
  # Get coverage probability over entire time axis (Overall)
  cov_prob_avg = mean(cov_prob_matrix)
  print(paste("Overall coverage probability:",cov_prob_avg))
  # Get coverage probability for each common time point
  cov_prob_over_time = rowMeans(cov_prob_matrix)
  
  # Convert to df to plot using ggplot
  cov_prob_over_time_df = data.frame(cov_prob_over_time)
  cov_prob_over_time_df$time = rownames(cov_prob_over_time_df)
  cov_prob_over_time_df$time = as.numeric(cov_prob_over_time_df$time)
  rownames(cov_prob_over_time_df) = NULL
  cov_prob_over_time_df$cov_prob_over_time = as.numeric(cov_prob_over_time_df$cov_prob_over_time)
  # print(paste("Coverage probability for each time point in common time axis:",cov_prob_over_time_df))
  
  cp_out_pref = "coverage_plots/"
  error_out_pref = "error_plots/"
  system(paste0("mkdir -p ",cp_out_pref,dirname(out_path_cov_prob)))
  system(paste0("mkdir -p ",error_out_pref,dirname(out_path_error)))
  system(paste0("mkdir -p ",error_out_pref,dirname(out_path_rmse)))

  p1 = ggplot(data=cov_prob_over_time_df, aes(x=time, y=cov_prob_over_time)) + 
    geom_line(color="steelblue") + geom_point(color="steelblue") + coord_cartesian(ylim = c(0,1)) + 
    scale_x_continuous(labels=scaleFUN, breaks = common_time_ax) + #scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Time axis", y = paste("Coverage probability for",length(sim_trees),"simulations"), caption = paste0("Overall coverage probability = ", scaleFUN(cov_prob_avg))) +
    theme_bw() + theme(plot.title = element_text(hjust=0.5), plot.caption = element_text(hjust=1),
                       axis.text.x = element_text(size=7, angle = 90, vjust = 0.5, hjust=1))

  print(paste0("Coverage probability for each time point plotted to: ",cp_out_pref,out_path_cov_prob))
  ggsave(paste0(cp_out_pref,out_path_cov_prob), plot=p1, units="in", width=7, height=7, dpi=600)
  
  p2 = ggplot(data=mean_abs_error_df, aes(x=as.numeric(as.character(n_sim)), y=mean_abs_error)) + 
    geom_line(color="steelblue") + geom_point(color="steelblue") +
    labs(x = "Simulation index", y = "Mean Absolute Error (MAE)", caption = paste0("Overall MAE = ", scaleFUN(mean_abs_error_avg))) +
    theme_bw() + theme(plot.title = element_text(hjust=0.5), plot.caption = element_text(hjust=1),
                       axis.text.x = element_text(size=7))
  
  print(paste0("Error (MAE) for each simulation index plotted to: ",error_out_pref,out_path_error))
  ggsave(paste0(error_out_pref,out_path_error), plot=p2, units="in", width=7, height=7, dpi=600)
  
  p3 = ggplot(data=rmse_df, aes(x=as.numeric(as.character(n_sim)), y=rmse_val)) + 
    geom_line(color="steelblue") + geom_point(color="steelblue") +
    labs(x = "Simulation index", y = "Root Mean Square Error (RMSE)", caption = paste0("Overall RMSE = ", scaleFUN(rmse_avg))) +
    theme_bw() + theme(plot.title = element_text(hjust=0.5), plot.caption = element_text(hjust=1),
                       axis.text.x = element_text(size=7))
  
  print(paste0("Root Mean Square Error (RMSE) for each simulation index plotted to: ",error_out_pref,out_path_rmse))
  ggsave(paste0(error_out_pref,out_path_rmse), plot=p3, units="in", width=7, height=7, dpi=600)
  
  return(list(pboot_ne_ind, pboot_ne_adj_df, cov_prob, cov_prob_matrix, cov_prob_avg, cov_prob_over_time, cov_prob_over_time_df, common_time_ax, mean_abs_error_df, rmse_df))
}

# Compare coverage probabilities for different models (skysigma, skygrid, and skygrowth) and the same Ne function
compare_cov_prob_models <- function(m1_cp, m2_cp, m3_cp, common_t_ax, out_path) {
  cp_out_pref = "coverage_plots/"
  system(paste0("mkdir -p ",cp_out_pref,dirname(out_path)))
  
  m1_cp$model = "Skysigma"; m2_cp$model = "Skygrid"; m3_cp$model = "Skygrowth"
  m_comb = rbind(m1_cp, m2_cp, m3_cp)
  m_comb$model = factor(m_comb$model, levels=c("Skysigma", "Skygrid", "Skygrowth"))
  
  p = ggplot(m_comb, aes(time, cov_prob_over_time, color = model)) + #shape=model
    geom_line() + geom_point() + scale_color_manual(name = "Model", values = c("steelblue", "darkred", "darkolivegreen")) +
    coord_cartesian(ylim = c(0,1)) +
    scale_x_continuous(labels=scaleFUN, breaks = common_t_ax) +
    labs(x = "Time axis", y = "Coverage probability for 500 simulations") +
    theme_bw() + theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(size=7, angle = 90, vjust = 0.5, hjust=1))
  
  print(paste0("Coverage probability for each time point plotted to: ",cp_out_pref,out_path))
  ggsave(paste0(cp_out_pref,out_path), plot=p, units="in", width=7, height=7, dpi=600)
  
  return(p)
}

# Compare coverage probability estimates for the same model using different sample sizes to see precision loss when dropping tips
compare_same_model_diff_samp_size <- function(cp_more_tips, cp_intermed_tips, cp_less_tips, common_t_ax, out_path) {
  cp_out_pref = "coverage_plots/"
  system(paste0("mkdir -p ",cp_out_pref,dirname(out_path)))
  
  cp_more_tips$sample_size = 200; cp_intermed_tips$sample_size = 50; cp_less_tips$sample_size = 20
  cp_comb = rbind(cp_more_tips, cp_intermed_tips, cp_less_tips)
  cp_comb$sample_size = factor(cp_comb$sample_size, levels=c(200, 50, 20))
  
  p = ggplot(cp_comb, aes(time, cov_prob_over_time, color = sample_size)) + 
    geom_line() + geom_point() + scale_color_manual(name = "Sample size", values = c("steelblue", "darkred", "darkolivegreen")) +
    coord_cartesian(ylim = c(0,1)) +
    scale_x_continuous(labels=scaleFUN, breaks = common_t_ax) +
    labs(x = "Time axis", y = "Coverage probability for 500 simulations") +
    theme_bw() + theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(size=7, angle = 90, vjust = 0.5, hjust=1))
  
  print(paste0("Coverage probability for each time point plotted to: ",cp_out_pref,out_path))
  ggsave(paste0(cp_out_pref,out_path), plot=p, units="in", width=7, height=7, dpi=600)
  
  return(p)
}

# Compare errors (MAE and RMSE) for different models (skysigma, skygrid, and skygrowth) and the same Ne function
compare_err_models <- function(m1_mae, m2_mae, m3_mae, m1_rmse, m2_rmse, m3_rmse, common_t_ax, out_path_mae, out_path_rmse) {
  err_out_pref = "error_plots/"
  system(paste0("mkdir -p ",err_out_pref,dirname(out_path_mae)))
  system(paste0("mkdir -p ",err_out_pref,dirname(out_path_rmse)))
  
  m1_mae$model = "Skysigma"; m2_mae$model = "Skygrid"; m3_mae$model = "Skygrowth"
  mae_comb = rbind(m1_mae, m2_mae, m3_mae)
  mae_comb$model = factor(mae_comb$model, levels=c("Skysigma", "Skygrid", "Skygrowth"))
  
  p1 = ggplot(mae_comb, aes(x=as.numeric(as.character(n_sim)), y=mean_abs_error, color = model)) +
    geom_line() + geom_point() + scale_color_manual(name = "Model", values = c("steelblue", "darkred", "darkolivegreen")) +
    labs(x = "Simulation index", y = "Mean Absolute Error (MAE)") +
    theme_bw() + theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(size=7))
  
  print(paste0("Mean absolute error (MAE) comparison for different models plotted to: ",err_out_pref,out_path_mae))
  ggsave(paste0(err_out_pref,out_path_mae), plot=p1, units="in", width=7, height=7, dpi=600)
  
  m1_rmse$model = "Skysigma"; m2_rmse$model = "Skygrid"; m3_rmse$model = "Skygrowth"
  rmse_comb = rbind(m1_rmse, m2_rmse, m3_rmse)
  rmse_comb$model = factor(rmse_comb$model, levels=c("Skysigma", "Skygrid", "Skygrowth"))
  
  p2 = ggplot(rmse_comb, aes(x=as.numeric(as.character(n_sim)), y=rmse_val, color = model)) +
    geom_line() + geom_point() + scale_color_manual(name = "Model", values = c("steelblue", "darkred", "darkolivegreen")) +
    labs(x = "Simulation index", y = "Root Mean Square Error (RMSE)") +
    theme_bw() + theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(size=7))
  
  print(paste0("Root Mean Square Error (RMSE) comparison for different models plotted to: ",err_out_pref,out_path_rmse))
  ggsave(paste0(err_out_pref,out_path_rmse), plot=p2, units="in", width=7, height=7, dpi=600)
  
  return(list(p1, p2))
}

# Compare error estimates (MAE and RMSE) for the same model using different sample sizes
compare_err_same_model_diff_samp_size <- function(mae_more_tips, mae_intermed_tips, mae_less_tips, rmse_more_tips, rmse_intermed_tips, rmse_less_tips, common_t_ax, out_path_mae, out_path_rmse) {
  err_out_pref = "error_plots/"
  system(paste0("mkdir -p ",err_out_pref,dirname(out_path_mae)))
  system(paste0("mkdir -p ",err_out_pref,dirname(out_path_rmse)))
  
  mae_more_tips$sample_size = 200; mae_intermed_tips$sample_size = 50; mae_less_tips$sample_size = 20
  mae_comb = rbind(mae_more_tips, mae_intermed_tips, mae_less_tips)
  mae_comb$sample_size = factor(mae_comb$sample_size, levels=c(200, 50, 20))
  
  p1 = ggplot(mae_comb, aes(x=as.numeric(as.character(n_sim)), y=mean_abs_error, color = sample_size)) + 
    geom_line() + geom_point() + scale_color_manual(name = "Sample size", values = c("steelblue", "darkred", "darkolivegreen")) +
    labs(x = "Simulation index", y = "Mean Absolute Error (MAE)") +
    theme_bw() + theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(size=7))
  
  print(paste0("Mean absolute error (MAE) comparison for same model and different sample sizes plotted to: ",err_out_pref,out_path_mae))
  ggsave(paste0(err_out_pref,out_path_mae), plot=p1, units="in", width=7, height=7, dpi=600)
  
  rmse_more_tips$sample_size = 200; rmse_intermed_tips$sample_size = 50; rmse_less_tips$sample_size = 20
  rmse_comb = rbind(rmse_more_tips, rmse_intermed_tips, rmse_less_tips)
  rmse_comb$sample_size = factor(rmse_comb$sample_size, levels=c(200, 50, 20))
  
  p2 = ggplot(rmse_comb, aes(x=as.numeric(as.character(n_sim)), y=rmse_val, color = sample_size)) + 
    geom_line() + geom_point() + scale_color_manual(name = "Sample size", values = c("steelblue", "darkred", "darkolivegreen")) +
    labs(x = "Simulation index", y = "Root Mean Square Error (RMSE)") +
    theme_bw() + theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(size=7))
  
  print(paste0("Root Mean Square Error (RMSE) comparison for same model and different sample sizes plotted to: ",err_out_pref,out_path_rmse))
  ggsave(paste0(err_out_pref,out_path_rmse), plot=p2, units="in", width=7, height=7, dpi=600)
  
  return(list(p1, p2))
}