#load task_table
load("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/nObsnGenTaskTable.rda")
n_tasks <- nrow(task_table)
n_evals <- 40
n_shapes <- length(unique(task_table$shape_name))

#make table to store results
sum_table <- data.frame("shape_name"= task_table$shape_name,
                        "nom_cover" = rep(c(.8, .9, .95), n_shapes),
                        "mean_cover_10" = rep(NA, 3*n_shapes),
                        "sd_cover_10"= rep(NA, 3*n_shapes),
                        "mean_cover_20" = rep(NA, 3*n_shapes),
                        "sd_cover_20"= rep(NA, 3*n_shapes),
                        "mean_cover_50" = rep(NA, 3*n_shapes),
                        "sd_cover_50"= rep(NA, 3*n_shapes))


for (i in 1:n_tasks) {
  n_obs_i <- task_table$n_obs[i]
  shape_name_i <- task_table$shape_name[i]
  n_gen_i <- task_table$n_gen[i]
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/nObsnGen/nObsnGen_id%i_%s_obs%i_gen%i.rda",
               i, shape_name_i, n_obs_i, n_gen_i))
  cover <- res_cover$cover/n_evals
  row_inds <- which(sum_table$shape_name == shape_name_i)
  if (n_obs_i == 10) {
    sum_table[row_inds, c("mean_cover_10")] <- apply(cover, 2, mean)
    sum_table[row_inds, c("sd_cover_10")] <- apply(cover, 2, sd)
  } else if (n_obs_i == 20) {
    sum_table[row_inds, c("mean_cover_20")] <- apply(cover, 2, mean)
    sum_table[row_inds, c("sd_cover_20")] <- apply(cover, 2, sd)
  } else {
    sum_table[row_inds, c("mean_cover_50")] <- apply(cover, 2, mean)
    sum_table[row_inds, c("sd_cover_50")] <- apply(cover, 2, sd)
  }
}

library("xtable")
xtable(sum_table, digits = 2)
