#Assess coverage under particular conditions
rm(list = ls())
set.seed(103)
library("ContouR")

#Read in simulation settings
task_id <- 3
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/variedPTaskTable.rda"
# task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/nObsnGenTaskTable.rda"
# task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/misspecTaskTable.rda"
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
                      sigma = sigma_true, C = C_true,thetas = thetas_true, 
                      r1_min = r1_min, r1_max = r1_max, r2_min = r2_min,
                      r2_max = r2_max, n_curl_min = n_curl_min,
                      n_curl_max = n_curl_max, bd = bd)
} else {
  test <- gen_conts(n_sim = n_evals, mu = mu_true, kappa = kappa_true,
                    sigma = sigma_true, C = C_true, thetas = thetas_true, 
                    bd = bd)
}

#grid of points for possible C
x_pts <- y_pts <- seq(0, 1, length = n_C_poss)
C_poss <- SpatialPoints(expand.grid(x_pts, y_pts))

#intial theta's with high p
theta_space_init <- 2*pi/p_init
thetas_init <- seq(theta_space_init/2, 2*pi, theta_space_init)

for (k in 1:n_evals) {
  start_time <- proc.time()  
  #simulate observations
  if (misspec) {
    obs <- gen_misspec(n_sim = n_obs, mu = mu_true, kappa = kappa_true,
                       sigma = sigma_true, C = C_true,thetas = thetas_true, 
                       r1_min = r1_min, r1_max = r1_max, r2_min = r2_min, 
                       r2_max = r2_max, n_curl_min = n_curl_min,
                       n_curl_max = n_curl_max, bd = bd)
  } else {
    obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true,
                     sigma = sigma_true, C = C_true, thetas = thetas_true,
                     bd = bd)
  }
  
  
  #find estimated center point, p, and theta
  area_tol <- err_prop*mean(sapply(obs$polys, gArea)) 
  C_est <- best_C(bd = bd, conts = obs$polys, thetas = thetas_init, 
                  area_tol = area_tol)
  p_est <- reduce_p(C = C_est, conts = obs$polys, area_tol = area_tol, p = p_init, 
                    red_prop = .05)
  theta_space_est <- 2*pi/p_est
  thetas_est <- seq(theta_space_est/2, 2*pi, theta_space_est)
  
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
  kappa_ini <- .1
  sigma_ini <-  apply(y, 1, sd)
  
  #priors
  mu0_indiv <- task$mu0
  mu0 <- rep(mu0_indiv, p_est)
  Lambda0_sigma2 <- task$Lambda0_sigma2
  Lambda0 <- Lambda0_sigma2*diag(p_est) 
  betaKappa0 <- task$betaKappa0
  betaSigma0 <- rep(task$betaSigma0, p_est)
  
  #MCMC set up
  g_space <- task$g_space
  g_space <- 1
  g_start <- seq(1, p_est, by = g_space)
  g_end <- c(seq(g_space, p_est, by = g_space), p_est)
  
  #non-scaler proposals for MCMC
  theta_dist_est <- theta_dist_mat(thetas_est)
  sigmaPropCov = sigmaProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)
  muPropCov <- muProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)
  
  #run MCMC
  fits <- ContouR::RunMCMC(nIter = n_iter, y,
                           mu = mu_ini, mu0 = mu0, Lambda0, muPropCov,
                           kappa = kappa_ini, betaKappa0, kappaProp_SD,
                           sigma = sigma_ini, betaSigma0 = betaSigma0, sigmaPropCov,
                           gStart = g_start - 1, gEnd = g_end - 1,
                           thetaDist = theta_dist_est)
  
  #parameter estimates
  mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
  kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
  sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)
  #rm(fits) #reduce memory
  
  #posterior field
  gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est,
                    sigma = sigma_est, C = C_est, thetas = thetas_est, bd = bd)
  prob <- prob_field(polys = gens$polys, nrows = n_grid, ncols = n_grid)
  #rm(gens) #reduce memory
  
  #find credible intervals and compute coverage
  creds <- cred_regs(prob, cred_eval = cred_levels, nrows = n_grid, 
                     ncols = n_grid)
  #rm(prob) #reduce memory
  if (k != 1) {
    cover <- cover + sapply(creds, 
                            function(x){eval_cred_reg(truth = test$polys[[k]],
                                                      cred_reg = x, 
                                                      center = C_true, 
                                                      p_test = p_test,
                                                      nrows = n_grid, ncols = n_grid)})
  } else {
    cover <- sapply(creds, 
                    function(x){eval_cred_reg(truth = test$polys[[k]],
                                              cred_reg = x, 
                                              center = C_true, 
                                              p_test = p_test,
                                              nrows = n_grid, ncols = n_grid)})
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
  file_name <- sprintf("%s_id%i_%s_pProp%f", task_name, task_id, shape_name, 
                       err_prop)
} else if (task_name == "misspec") {
  file_name <- sprintf("%s_id%i_%s_nCurlMin%i", task_name, task_id, shape_name,
                       n_curl_min)
}