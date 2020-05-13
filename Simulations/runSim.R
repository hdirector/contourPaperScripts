#Assess coverage under particular conditions
rm(list = ls())
set.seed(101)
library("ContouR")

#Read in arguments
args <- commandArgs(trailingOnly =TRUE)
task_id <- as.numeric(args[1])
task_name <- args[2]
task_path <- sprintf("/homes/direch/contours/Simulations/%sTaskTable.rda", task_name)
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
                      n_curl_max = n_curl_max, bd = bd,  rand_loc = rand_loc)
} else {
  test <- gen_conts(n_sim = n_evals, mu = mu_true, kappa = kappa_true,
                    sigma = sigma_true, C = C_true, thetas = thetas_true, 
                    bd = bd)
}

#storage for p_est and C_est
p_est_store <- rep(NA, n_evals)
C_est_store <- matrix(nrow = n_evals, ncol = 2, data = NA)

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
  if (!misspec) {
	ests <- find_CP(conts = obs, delta, p_init, space, step, misspec)
	C_est <- ests$C_est
	thetas_est <- ests$thetas_est
	p_est <- length(thetas_est)
  } else {
    p_est <- p_init
    theta_space_est <- 2*pi/p_est
    thetas_est <- seq(theta_space_est/2, 2*pi, theta_space_est)
    C_est <- best_C(bd = bd, conts = obs$polys, thetas = thetas_est,
                    space = space)
  }
  C_est_store[k,] <- C_est
  p_est_store[k] <- p_est

  #measure and store y
  #Make sets of lines, l, for  C_hat and get y's
  l_untrim <- make_l(C = C_est, thetas_est)
  l <- lapply(l_untrim, function(x){gIntersection(x, box)})
  l_lengths <- sapply(l, function(x){as.numeric(gLength(x))})
  
  #compute y's
  pts <- lapply(obs$polys, function(x) {pts_on_l(l, x, under = FALSE)})
  y <- sapply(pts, function(x){apply(x, 1, function(y){get_dist(y, C_est)})})
  rm(obs) #reduce memory
  
  #initial values for MCMC
  mu_ini <- apply(y, 1, mean)
  kappa_ini <- 1
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
  rm(fits) #reduce memory
 
  #posterior field
  gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est,
                    sigma = sigma_est, C = C_est, thetas = thetas_est, bd = bd)
  prob <- prob_field(polys = gens$polys, nrows = n_grid, ncols = n_grid)
  rm(gens) #reduce memory
  
  #find credible intervals and compute coverage
  creds <- cred_regs(prob, cred_eval = cred_levels, nrows = n_grid, 
                     ncols = n_grid)
  
  #compute theta's to evaluate
  theta_space_eval <- 2*pi/p_eval
  thetas_eval <- seq(theta_space_eval/2, 2*pi, by = theta_space_eval)
  
  if (k != 1) {
    cover <- cover + sapply(creds, 
                            function(x){eval_cred_reg(truth = test$polys[[k]],
                                                      cred_reg = x, 
                                                      center = C_true, 
                                                      thetas = thetas_eval,
                                                      nrows = n_grid, 
                                                      ncols = n_grid)})
  } else {
    cover <- sapply(creds, 
                    function(x){eval_cred_reg(truth = test$polys[[k]],
                                              cred_reg = x, 
                                              center = C_true, 
                                              thetas = thetas_eval,
                                              nrows = n_grid, ncols = n_grid)})
  }
  rm(creds) #reduce memory
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
} else if (task_name == "shapes") {
  file_name <- sprintf("%s_id%i_%s", task_name, task_id, shape_name)
}
res <- list("res_cover" = res_cover, "p_est" = p_est_store, "C_est" = C_est_store)
save(res, file = sprintf("/homes/direch/contours/Simulations/sim_results/%s/%s.rda", task_name, file_name))
