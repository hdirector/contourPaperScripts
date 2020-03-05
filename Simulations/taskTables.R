#########################################################################
#Simulate varied number of observations and generated contours
#########################################################################
#varied settings
n_obs <- c(5, 20, 50)
n_gen <- c(25, 100, 200)
shape_name <- c("shape2_2")
task_table <- as.data.frame(expand.grid(n_obs, n_gen, shape_name))
colnames(task_table) <- c("n_obs", "n_gen", "shape_name")
n_case <- nrow(task_table)

#fixed settings
task_table$p_test <- rep(100, n_case)
task_table$p_true <- rep(50, n_case)
task_table$p_init <- rep(50, n_case)
task_table$n_evals <- rep(50, n_case) 
task_table$n_grid <- rep(200, n_case)
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
task_table$err_prop <- rep(.025, n_case)
task_table$misspec <- rep(FALSE, n_case)
task_table$r1_min <- rep(NULL, n_case)
task_table$r1_max <- rep(NULL, n_case)
task_table$r2_min <- rep(NULL, n_case)
task_table$r2_max <- rep(NULL, n_case)
task_table$n_curl_min <- rep(NULL, n_case)
task_table$n_curl_max <- rep(NULL, n_case)
save(task_table, 
     file = '/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/nObsnGenTaskTable.rda')


#######################################################################
#Simulate varied proportions of the true p that are used 
#######################################################################
#varied settings
err_prop <- c(.05, .025, .01)
shape_name <- c("shape2_1", "shape2_2", "shape2_5")
task_table <- as.data.frame(expand.grid(err_prop, shape_name))
colnames(task_table) <- c("err_prop",  "shape_name")
n_case <- nrow(task_table)

#fixed settings
task_table$n_obs <- rep(20, n_case)
task_table$n_gen <- rep(100, n_case)
task_table$p_test <- rep(100, n_case)
task_table$p_true <- rep(50, n_case)
task_table$p_init <- rep(50, n_case)
task_table$n_evals <- rep(50, n_case) 
task_table$n_grid <- rep(200, n_case) 
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
task_table$misspec <- rep(FALSE, n_case)
task_table$r1_min <- rep(NULL, n_case)
task_table$r1_max <- rep(NULL, n_case)
task_table$r2_min <- rep(NULL, n_case)
task_table$r2_max <- rep(NULL, n_case)
task_table$n_curl_min <- rep(NULL, n_case)
task_table$n_curl_max <- rep(NULL, n_case)

save(task_table, 
     file = '/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/variedPTaskTable.rda')



#############################################
#Simulate fits under model misspecification
#############################################
#varied settings
shape_name <- c("shape2_2")
n_curl_min <- c(1, 2, 3)
task_table <- as.data.frame(expand.grid(n_curl_min, shape_name))
colnames(task_table) <- c("n_curl_min", "shape_name")
n_case <- nrow(task_table)

#fixed settings
task_table$n_obs <- rep(20, n_case)
task_table$n_gen <- rep(100, n_case)
task_table$p_true <- rep(50, n_case)
task_table$p_test <- rep(100, n_case)
task_table$n_evals <- rep(50, n_case)
task_table$n_grid <- rep(200, n_case)
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
task_table$err_prop <- rep(.1, n_case)
task_table$misspec <- rep(TRUE, n_case)
task_table$r1_min <- rep(0.02, n_case)
task_table$r1_max <- rep(.04, n_case)
task_table$r2_min <- rep(.05, n_case)
task_table$r2_max <- rep(.1, n_case)
task_table$n_curl_max <- task_table$n_curl_min + 1
save(task_table,
     file = '/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/misspecTaskTable.rda')

