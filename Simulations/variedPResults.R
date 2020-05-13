rm(list = ls())
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/variedPTaskTable.rda"

load(task_path)
n_tasks <- nrow(task_table)
n_eval <- 40

res_variedP <- data.frame(cbind(rep(c(1, 2, 4), each = 3), rep(c(.8, .9, .95),3),
                     matrix(nrow = 9, ncol = 12)))
colnames(res_variedP) <- c("kappa", "nom_cov", "mean_03", "sd_03", "p_03", "sdp_03", 
                           "mean_02", "sd_02", "p_02", "sdp_02",  
                           "mean_01", "sd_01", "p_01", "sdp_01")

for (i in 1:n_tasks) {
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/variedP/variedP_id%i_%s_delta%f.rda",
               i, task_table[i,]$shape, task_table[i,]$delta))
  if (task_table$shape_name[i] == "shapeA_1") {
    row_ind <- 1:3 
  } else if (task_table$shape_name[i] == "shapeA_2") {
    row_ind <- 4:6
  } else {
    row_ind <- 7:9
  }
  if (task_table$delta[i] == .03) {
    res_variedP[row_ind, "mean_03"] <-  apply(res$res_cover$cover, 2, mean)/n_eval
    res_variedP[row_ind, "sd_03"] <-  apply(res$res_cover$cover/n_eval, 2, sd)
    res_variedP[row_ind, "p_03"] <- rep(mean(res$p_est), 3)
    res_variedP[row_ind, "sdp_03"] <- rep(sd(res$p_est), 3)
  } else if (task_table$delta[i] == .02) {
    res_variedP[row_ind, "mean_02"] <-  apply(res$res_cover$cover, 2, mean)/n_eval
    res_variedP[row_ind, "sd_02"] <-  apply(res$res_cover$cover/n_eval, 2, sd)
    res_variedP[row_ind, "p_02"] <- rep(mean(res$p_est), 3)
    res_variedP[row_ind, "sdp_02"] <- rep(sd(res$p_est), 3)
  } else {
    res_variedP[row_ind, "mean_01"] <-  apply(res$res_cover$cover, 2, mean)/n_eval
    res_variedP[row_ind, "sd_01"] <-  apply(res$res_cover$cover/n_eval, 2, sd)
    res_variedP[row_ind, "p_01"] <- rep(mean(res$p_est), 3)
    res_variedP[row_ind, "sdp_01"] <- rep(sd(res$p_est), 3)
  }
}

rownames(res_variedP) <- c("0.8", "0.9", "0.95")

xtable::xtable(res_variedP, digits = 2)

