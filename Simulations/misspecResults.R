rm(list = ls())
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/misspecTaskTable.rda"
load(task_path)

n_tasks <- nrow(task_table)
n_evals <- 50

res_misspec <- data.frame("mean_0_1" = rep(NA, 3), "sd_0_1" = rep(NA, 3),
                         "mean_2_3" = rep(NA, 3), "sd_2_3" = rep(NA, 3),
                         "mean_4_5" = rep(NA, 3),"sd_4_5" = rep(NA, 3))

rownames(res_misspec)[1:3] <- c("0.8", "0.9", "0.95")


for (i in 1:n_tasks) {
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/misspec/misspec_id%i_%s_nCurlMin%i.rda",
               i, task_table[i,]$shape, task_table[i,]$n_curl_min))
  if (task_table[i,]$n_curl_min == 0) {
    res_misspec[,"mean_0_1"] <- apply(res_cover$cover, 2, mean)/n_evals
    res_misspec[,"sd_0_1"] <- apply(res_cover$cover, 2, sd)/n_evals
  } else if (task_table[i,]$n_curl_min == 2) {
    res_misspec[,"mean_2_3"] <- apply(res_cover$cover, 2, mean)/n_evals
    res_misspec[,"sd_2_3"] <- apply(res_cover$cover, 2, sd)/n_evals
  } else if (task_table[i,]$n_curl_min == 4) {
    res_misspec[,"mean_4_5"] <- apply(res_cover$cover, 2, mean)/n_evals
    res_misspec[,"sd_4_5"] <- apply(res_cover$cover, 2, sd)/n_evals
  }
    
}

xtable::xtable(res_misspec, digits = 3)
