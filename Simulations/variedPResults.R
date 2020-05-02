rm(list = ls())
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/variedPTaskTable.rda"

load(task_path)
n_tasks <- nrow(task_table)
n_eval <- 50

res_variedP <- data.frame("mean_05" = rep(NA, 3), "sd_05" = rep(NA, 3),
                          "mean_03" = rep(NA, 3), "sd_03" = rep(NA, 3),
                          "mean_01" = rep(NA, 3),"sd_01" = rep(NA, 3))

for (i in 1:n_tasks) {
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/variedP/variedP_id%i_%s_delta%f.rda",
               i, task_table[i,]$shape, task_table[i,]$delta))
  if (task_table[i,]$delta == .05) {
    res_variedP$mean_05 <-  apply(res_cover$cover, 2, mean)/n_eval
    res_variedP$sd_05 <- apply(res_cover$cover/n_eval, 2, sd)
  } else if (task_table[i,]$delta == .03) {
    res_variedP$mean_03 <-  apply(res_cover$cover, 2, mean)/n_eval
    res_variedP$sd_03 <- apply(res_cover$cover/n_eval, 2, sd)
  } else if (task_table[i,]$delta == .01) {
    res_variedP$mean_01 <-  apply(res_cover$cover, 2, mean)/n_eval
    res_variedP$sd_01 <- apply(res_cover$cover/n_eval, 2, sd)
  }
}

rownames(res_variedP) <- c("0.8", "0.9", "0.95")

xtable::xtable(res_variedP, digits = 3)

