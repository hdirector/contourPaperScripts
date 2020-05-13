rm(list = ls())
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/shapesTaskTable.rda"

load(task_path)
n_tasks <- nrow(task_table)
n_eval <- 40

res_shapes <- data.frame("mean_A" = rep(NA, 3), "sd_A" = rep(NA, 3),
                          "mean_B" = rep(NA, 3), "sd_B" = rep(NA, 3),
                          "mean_C" = rep(NA, 3),"sd_C" = rep(NA, 3))

for (i in 1:n_tasks) {
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/shapes/shapes_id%i_%s.rda",
               i, task_table[i,]$shape))
  if (task_table[i,]$shape_name == "shapeA_2") {
    res_shapes$mean_A <-  apply(res_cover$cover, 2, mean)/n_eval
    res_shapes$sd_A <- apply(res_cover$cover/n_eval, 2, sd)
  } else if (task_table[i,]$shape_name == "shapeB") {
    res_shapes$mean_B <-  apply(res_cover$cover, 2, mean)/n_eval
    res_shapes$sd_B <- apply(res_cover$cover/n_eval, 2, sd)
  } else if (task_table[i,]$shape_name == "shapeC") {
    res_shapes$mean_C <-  apply(res_cover$cover, 2, mean)/n_eval
    res_shapes$sd_C <- apply(res_cover$cover/n_eval, 2, sd)
  }
}

rownames(res_shapes) <- c("0.8", "0.9", "0.95")

library("xtable")
xtable(res_shapes, digits = 2)
