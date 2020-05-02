#########################################################################
#Simulate varied num obs, num gen, and correlation kappa
#########################################################################
#varied settings
n_obs <- c(10, 20, 50, 100)
n_gen <- 100
shape_name <- c("shapeA_1", "shapeA_2", "shapeA_4")
task_table <- as.data.frame(expand.grid(n_obs, n_gen, shape_name))
colnames(task_table) <- c("n_obs", "n_gen", "shape_name")
n_case <- nrow(task_table)

#fixed settings
task_table$under <- rep(TRUE, n_case)
task_table$p_eval <- rep(100, n_case)
task_table$p_true <- rep(50, n_case)
task_table$p_init <- rep(46, n_case)
task_table$n_evals <- rep(50, n_case) 
task_table$n_grid <- rep(200, n_case)
task_table$mu0 <- rep(.2, n_case)
task_table$Lambda0_sigma2 <- rep(.05, n_case)
task_table$betaKappa0 <- rep(8, n_case)
task_table$betaSigma0 <- rep(.15, n_case)
task_table$n_iter <- rep(50000, n_case)
task_table$burn_in <- rep(15000, n_case)
task_table$sigmaProp_sigma2 <- rep(.05, n_case)
task_table$muProp_sigma2 <- rep(.05, n_case)
task_table$kappaProp_SD <- rep(.3, n_case)
task_table$delta <- rep(.01, n_case)
task_table$space <- rep(.01, n_case)
task_table$step <- rep(.05, n_case)
task_table$misspec <- rep(FALSE, n_case)
task_table$r1_min <- rep(NULL, n_case)
task_table$r1_max <- rep(NULL, n_case)
task_table$r2_min <- rep(NULL, n_case)
task_table$r2_max <- rep(NULL, n_case)
task_table$n_curl_min <- rep(NULL, n_case)
task_table$n_curl_max <- rep(NULL, n_case)
task_table$rand_loc <- rep(NULL, n_case)
save(task_table, 
     file = '/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/nObsnGenTaskTable.rda')


######################################################
#Use different delta's 
######################################################
#varied settings
delta <- c(.05, .025, .01)
shape_name <- c("shapeA_2")
task_table <- as.data.frame(expand.grid(delta, shape_name))
colnames(task_table) <- c("delta",  "shape_name")
n_case <- nrow(task_table)

#fixed settings
task_table$under <- rep(TRUE, n_case)
task_table$n_obs <- rep(20, n_case)
task_table$n_gen <- rep(100, n_case)
task_table$p_eval <- rep(100, n_case)
task_table$p_true <- rep(50, n_case)
task_table$p_init <- rep(30, n_case)
task_table$n_evals <- rep(50, n_case) 
task_table$n_grid <- rep(200, n_case) 
task_table$mu0 <- rep(.2, n_case)
task_table$Lambda0_sigma2 <- rep(.05, n_case)
task_table$betaKappa0 <- rep(8, n_case)
task_table$betaSigma0 <- rep(.15, n_case)
task_table$n_iter <- rep(50000, n_case)
task_table$burn_in <- rep(15000, n_case)
task_table$sigmaProp_sigma2 <- rep(.05, n_case)
task_table$muProp_sigma2 <- rep(.05, n_case)
task_table$kappaProp_SD <- rep(.3, n_case)
task_table$space <- rep(.01, n_case)
task_table$step <- rep(.05, n_case)
task_table$misspec <- rep(FALSE, n_case)
task_table$r1_min <- rep(NULL, n_case)
task_table$r1_max <- rep(NULL, n_case)
task_table$r2_min <- rep(NULL, n_case)
task_table$r2_max <- rep(NULL, n_case)
task_table$n_curl_min <- rep(NULL, n_case)
task_table$n_curl_max <- rep(NULL, n_case)
task_table$rand_loc <- rep(NULL, n_case)


save(task_table, 
     file = '/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/variedPTaskTable.rda')

#############################################
#Simulate fits under model misspecification
#############################################
#varied settings
shape_name <- c("shapeA_2")
n_curl_min <- c(0, 2, 4)
rand_loc <- c("TRUE", "FALSE")
task_table <- as.data.frame(expand.grid(n_curl_min, shape_name,
                                        rand_loc))
colnames(task_table) <- c("n_curl_min", "shape_name", "rand_loc")
n_case <- nrow(task_table)

#fixed settings
task_table$under <- rep(TRUE, n_case)
task_table$n_obs <- rep(20, n_case)
task_table$n_gen <- rep(100, n_case)
task_table$p_true <- rep(50, n_case)
task_table$p_eval <- rep(100, n_case)
task_table$p_init <- rep(46, n_case)
task_table$n_evals <- rep(50, n_case)
task_table$n_grid <- rep(200, n_case)
task_table$mu0 <- rep(.2, n_case)
task_table$Lambda0_sigma2 <- rep(.05, n_case)
task_table$betaKappa0 <- rep(8, n_case)
task_table$betaSigma0 <- rep(.15, n_case)
task_table$n_iter <- rep(50000, n_case)
task_table$burn_in <- rep(15000, n_case)
task_table$sigmaProp_sigma2 <- rep(.05, n_case)
task_table$muProp_sigma2 <- rep(.05, n_case)
task_table$kappaProp_SD <- rep(.3, n_case)
task_table$delta <- rep(.03, n_case)
task_table$misspec <- rep(TRUE, n_case)
task_table$space <- rep(.01, n_case)
task_table$step <- rep(.05, n_case)
task_table$r1_min <- rep(0.05, n_case)
task_table$r1_max <- rep(.075, n_case)
task_table$r2_min <- rep(.075, n_case)
task_table$r2_max <- rep(.15, n_case)
task_table$n_curl_max <- task_table$n_curl_min + 1
save(task_table,
     file = '/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/misspecTaskTable.rda')

#################################
#look at different shapes
#################################
#varied settings
task_table <- data.frame(shape_name =  c("shapeA_2", "shapeB", "shapeC"))
n_case <- nrow(task_table)

#fixed settings
task_table$under <- rep(TRUE, n_case)
task_table$n_obs <- rep(20, n_case)
task_table$n_gen <- rep(100, n_case)
task_table$p_true <- rep(50, n_case)
task_table$p_eval <- rep(100, n_case)
task_table$p_init <- rep(46, n_case)
task_table$n_evals <- rep(50, n_case)
task_table$n_grid <- rep(200, n_case)
task_table$mu0 <- rep(.2, n_case)
task_table$Lambda0_sigma2 <- rep(.05, n_case)
task_table$betaKappa0 <- rep(8, n_case)
task_table$betaSigma0 <- rep(.15, n_case)
task_table$n_iter <- rep(50000, n_case)
task_table$burn_in <- rep(15000, n_case)
task_table$sigmaProp_sigma2 <- rep(.05, n_case)
task_table$muProp_sigma2 <- rep(.05, n_case)
task_table$kappaProp_SD <- rep(.3, n_case)
task_table$delta <- rep(.01, n_case)
task_table$space <- rep(.01, n_case)
task_table$step <- rep(.05, n_case)
task_table$misspec <- rep(FALSE, n_case)
task_table$r1_min <- rep(NULL, n_case)
task_table$r1_max <- rep(NULL, n_case)
task_table$r2_min <- rep(NULL, n_case)
task_table$r2_max <- rep(NULL, n_case)
task_table$n_curl_max <- rep(NULL, n_case)
task_table$rand_loc <- rep(NULL, n_case)
save(task_table,
     file = '/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/shapesTaskTable.rda')


