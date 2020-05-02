#Assess coverage under particular conditions
rm(list = ls())
set.seed(103)
library("ContouR")
library("coda")

#Read in simulation settings
task_id <- 3
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/shapesTaskTable.rda"
#task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/variedPTaskTable.rda"
#task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/nObsnGenTaskTable.rda"
#task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/misspecTaskTable.rda"
load(task_path)
task <- task_table[task_id,]
print(sprintf("task_id: %i", task_id))
attach(task)

#credible intervals
cred_levels <- c(80, 90, 95)
n_cred <- length(cred_levels)

#true parameters
theta_space_true <- 2*pi/p_true
thetas_true <- seq(theta_space_true/2, 2*pi, theta_space_true)
shape_name <- as.character(task$shape_name)
pars <- get(shape_name)
C_true <- pars$C
mu_true <- pars$mu
sigma_true <- pars$sigma
kappa_true <- pars$kappa

#sample test contours from the truth for later testing
box <- bbox()
bd <- box@polygons[[1]]@Polygons[[1]]@coords
if (misspec) {
  test <- gen_misspec(n_sim = n_evals, mu = mu_true, kappa = kappa_true,
                      sigma = sigma_true, C = C_true, thetas = thetas_true, 
                      r1_min = r1_min, r1_max = r1_max, r2_min = r2_min,
                      r2_max = r2_max, n_curl_min = n_curl_min,
                      n_curl_max = n_curl_max, bd = bd)
} else {
  test <- gen_conts(n_sim = n_evals, mu = mu_true, kappa = kappa_true,
                    sigma = sigma_true, C = C_true, thetas = thetas_true, 
                    bd = bd)
}

for (k in 1:n_evals) {
  start_time <- proc.time()  
  #simulate observations
  if (misspec) {
    obs <- gen_misspec(n_sim = n_obs, mu = mu_true, kappa = kappa_true,
                       sigma = sigma_true, C = C_true,thetas = thetas_true, 
                       r1_min = r1_min, r1_max = r1_max, r2_min = r2_min, 
                       r2_max = r2_max, n_curl_min = n_curl_min,
                       n_curl_max = n_curl_max, bd = bd, rand_loc = rand_loc)
  } else {
    obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true,
                     sigma = sigma_true, C = C_true, thetas = thetas_true,
                     bd = bd)
  }
  
  
  #find estimated center point, p, and theta
  ests <- find_CP(conts = obs, delta, p_init, space, step, misspec,
                  parallel = TRUE)
  C_est <- ests$C_est
  thetas_est <- ests$thetas_est
  p_est <- ests$p_est
  
  #measure and store y
  #Make sets of lines, l, for  C_hat and get y's
  l_untrim <- make_l(C = C_est, thetas_est)
  l <- lapply(l_untrim, function(x){gIntersection(x, box)})
  l_lengths <- sapply(l, function(x){as.numeric(gLength(x))})
  
  #compute y's
  pts <- lapply(obs$polys, function(x) {pts_on_l(l, x, under = FALSE)})
  y <- sapply(pts, function(x){apply(x, 1, function(y){get_dist(y, C_est)})})
  #rm(obs) #reduce memory
  
  #initial values for MCMC
  mu_ini <- apply(y, 1, mean)
  kappa_ini <- 3
  sigma_ini <-  apply(y, 1, sd)
  
  #priors
  mu0_indiv <- task$mu0
  mu0 <- rep(mu0_indiv, p_est)
  Lambda0 <- Lambda0_sigma2*diag(p_est) 
  betaSigma0 <- rep(task$betaSigma0, p_est)
  
  #non-scaler proposals for MCMC
  theta_dist_est <- theta_dist_mat(thetas_est)
  sigmaPropSD = sqrt(diag(sigmaProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)))
  muPropSD <- sqrt(diag(muProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)))
  
  #run MCMC
  start_time <- proc.time()
  fits <- RunMCMC(nIter = n_iter, y,
                  mu = mu_ini, mu0 = mu0, Lambda0, muPropSD = muPropSD,
                  kappa = kappa_ini, betaKappa0, kappaProp_SD,
                  sigma = sigma_ini, betaSigma0 = betaSigma0, 
                  sigmaPropSD = sigmaPropSD,
                  thetaDist = theta_dist_est)
  end_time <- proc.time()
  run_time <- end_time - start_time
  print(sprintf("fitted eval %i, time: %s", k, run_time[3]))
  
  #parameter estimates
  mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
  kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
  sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)
  #rm(fits) #reduce memory
  
  par(mfrow = c(1, 3))
  Sigma_obs <- cov(t(y))
  Sigma_ini <- compSigma(sigma_ini, kappa_ini, theta_dist_est)
  Sigma_est <- compSigma(sigma_est, kappa_est, theta_dist_est)
  Sigma_true <- compSigma(sigma_true, kappa_true, theta_dist_est)
  image.plot(Sigma_est, zlim  = c(0, .007))
  image.plot(Sigma_obs, zlim = c(0, .007))
  image.plot(Sigma_true, zlim = c(0, .007))

  par(mfrow = c(1, 1))
  plot(sigma_est, type= "l")
  points(sigma_ini, type= "l", col = "blue")
  points(sigma_true, type= "l", col = 'red')
  points(betaSigma0, type = "l", lty = 2)
  plot(mu_est, type= "l")
  points(mu_ini, col= 'blue', type= "l")
  points(mu_true, type= 'l',col = 'red')

  # #check MCMC diagnostic
  # N_sigma <- N_mu <- rep(NA, p_est)
  # for (i in 1:p_est) {
  #   N_sigma[i] <- max(raftery.diag(fits$sigma[i,], q = .025, r = .0125)$resmatrix[,"N"],
  #                     raftery.diag(fits$sigma[i,], q = .975, r = .0125)$resmatrix[,"N"])
  #   N_mu[i] <- max(raftery.diag(fits$mu[i,], q = .025, r = .0125)$resmatrix[,"N"],
  #                   raftery.diag(fits$mu[i,], q = .975, r = .0125)$resmatrix[,"N"])
  # }
  # N_kappa <- max(raftery.diag(fits$kappa, q = .025, r = .0125)$resmatrix[,"N"],
  #                raftery.diag(fits$kappa, q = .975, r = .0125)$resmatrix[,"N"])

  #posterior field
  gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est,
                    sigma = sigma_est, C = C_est, thetas = thetas_est, bd = bd)
  prob <- prob_field(polys = gens$polys, nrows = n_grid, ncols = n_grid)
  #rm(gens) #reduce memory
  
  #find credible intervals and compute coverage
  creds <- cred_regs(prob, cred_eval = cred_levels, nrows = n_grid, 
                     ncols = n_grid)
  
  #compute theta's to evaluate
  theta_space_eval <- 2*pi/p_eval
  thetas_eval <- seq(theta_space_eval/2, 2*pi, by = theta_space_eval)
  
  par(mfrow = c(1, 3))
  if (k != 1) {
    cover <- cover + sapply(creds, 
                            function(x){eval_cred_reg(truth = test$polys[[k]],
                                                      cred_reg = x, 
                                                      center = C_true, 
                                                      thetas = thetas_eval,
                                                      nrows = n_grid, 
                                                      ncols = n_grid,
                                                      plotting = TRUE)})
    k <- k + 1
  } else {
    cover <- sapply(creds, 
                    function(x){eval_cred_reg(truth = test$polys[[k]],
                                              cred_reg = x, 
                                              center = C_true, 
                                              thetas = thetas_eval,
                                              nrows = n_grid, ncols = n_grid,
                                              plotting = TRUE)})
  }
  #rm(creds) #reduce memory
  end_time <- proc.time()
  elapse_time <- end_time  - start_time
  print(sprintf("Eval %i completed for task_id %i", k, task_id))
  print(elapse_time)
}

#save results
res_cover <- list("task" = task, "cover" = cover)
if (task_name == "nObsnGen") {
  file_name <- sprintf("%s_id%i_%s_obs%i_gen%i", task_name, task_id, shape_name,
                       n_obs, n_gen)
} else if (task_name == "variedP") {
  file_name <- sprintf("%s_id%i_%s_delta%f", task_name, task_id, shape_name, 
                       delta)
} else if (task_name == "misspec") {
  file_name <- sprintf("%s_id%i_%s_nCurlMin%i", task_name, task_id, shape_name,
                       n_curl_min)
}

