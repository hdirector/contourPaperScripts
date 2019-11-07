################################################################
#Simulate varied number of observations and generated contours
###############################################################
#varied settings
n_obs <- c(5, 15, 25, 50)
n_gen <- c(25, 50, 100, 200)
shape_name <- c("stop_sign", "tie", "tree")
task_table <- as.data.frame(expand.grid(n_obs, n_gen, shape_names))
colnames(task_table) <- c("n_obs", "n_gen", "shape_name")
n_case <- nrow(task_table)

#fixed settings
task_table$p_test <- rep(30, n_case)
task_table$n_evals <- rep(100, n_case) 
task_table$n_grid <- rep(300, n_case)
task_table$mu0 <- rep(.1, n_case)
task_table$Lambda0_sigma2 <- rep(.05, n_case)
task_table$betaKappa0 <- rep(100, n_case)
task_table$betaSigma0 <- rep(.1, n_case)
task_table$n_iter <- rep(30000, n_case)
task_table$burn_in <- rep(20000, n_case)
task_table$sigmaProp_sigma2 <- rep(.05, n_case)
task_table$muProp_sigma2 <- rep(.05, n_case)
task_table$kappaProp_SD <- rep(.3, n_case)
task_table$g_space <- rep(1, n_case)


