rm(list = ls())
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/misspecTaskTable.rda"
load(task_path)

n_tasks <- nrow(task_table)

res_misspec <- cbind(task_table[1:n_tasks,c("shape_name", "n_curl_min", "n_curl_max")],
             rep(NA, n_tasks), rep(NA, n_tasks), rep(NA, n_tasks))
colnames(res_misspec)[4:6] <- c("Nominal .8", "Nominal .9", "Nominal .95")


for (i in c(1:2, 4:n_tasks)) {
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/misspec/misspec_id%i_%s_nCurlMin%i.rda",
               i, task_table[i,]$shape, task_table[i,]$n_curl_min))
  res_misspec[i,  c("Nominal .8", "Nominal .9", "Nominal .95")] <- apply(res_cover$cover, 2, mean)/50
}

xtable::xtable(res_misspec, digits = 3)
