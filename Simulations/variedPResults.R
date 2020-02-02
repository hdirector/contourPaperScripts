rm(list = ls())
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/variedPTaskTableOld.rda"
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/variedPTaskTable.rda"

load(task_path)

res_variedP <- cbind(task_table[1:12,c("shape_name", "p_prop")], rep(NA, 12), rep(NA, 12), rep(NA, 12))
colnames(res_variedP)[3:5] <- c("Nominal .8", "Nominal .9", "Nominal .95")


for (i in 1:nrow(task_table)) {
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/variedP/variedP_id%i_%s_pProp%f.rda",
               i, task_table[i,]$shape, task_table[i,]$p_prop))
  res_variedP[i,  c("Nominal .8", "Nominal .9", "Nominal .95")] <- apply(res_cover$cover, 2, mean)/50
}

xtable::xtable(res_variedP, digits = 3)
