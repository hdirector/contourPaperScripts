#Assess coverage under particular conditions
rm(list = ls())
set.seed(103)
library("ContouR")

#Read in arguments
#args <- commandArgs(trailingOnly =TRUE)
task_id <- 19# as.numeric(args[1])
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/nObsnGenTaskTable.rda"
load(task_path)
task <- task_table[task_id,]
print(sprintf("task_id: %i", task_id))

#simulation settings
n_obs <- task$n_obs
n_gen <- task$n_gen
cred_levels <- c(80, 90, 95)
n_cred <- length(cred_levels)
p_test <- task$p_test
n_evals <- 10 #task$n_evals 
n_grid <- 100 #task$n_grid

#true parameters
pars <- get(as.character(task$shape_name))
p <- length(pars$mu)
mu_true <- pars$mu
p <- length(mu_true)
kappa_true <- .1 # pars$kappa
sigma_true <- pars$sigma
Cx_true <- pars$Cx
Cy_true <- pars$Cy
thetas_true <- pars$theta

#priors
mu0_indiv <- .3#task$mu0
mu0 <- rep(mu0_indiv, p)
Lambda0_sigma2 <- .05#task$Lambda0_sigma2
Lambda0 <- Lambda0_sigma2*diag(p) 
betaKappa0 <- 200 #task$betaKappa0
betaSigma0 <- rep(task$betaSigma0, p)

#MCMC settings
n_iter <- task$n_iter
burn_in <- task$burn_in
sigmaProp_sigma2 <- task$sigmaProp_sigma2
muProp_sigma2 <- task$muProp_sigma2
kappaPropSD <- task$kappaProp_SD
g_space <- task$g_space
g_space <- 1
g_start <- seq(1, p, by = g_space)
g_end <- c(seq(g_space, p, by = g_space), p)

#sample test contours from the truth for later testing
box <- bbox()
bd <- box@polygons[[1]]@Polygons[[1]]@coords
test <- gen_conts(n_sim = n_evals, mu = mu_true, kappa = kappa_true,
                  sigma = sigma_true, Cx = Cx_true, Cy = Cy_true,
                  thetas = thetas_true, bd = bd)

#grid of points for possible C
n_C_poss <- 20
x_pts <- y_pts <- seq(0, 1, length = n_C_poss)
C_poss <- SpatialPoints(expand.grid(x_pts, y_pts))


#set up thetas
theta_space <- 2*pi/p
thetas <- seq(theta_space, 2*pi, theta_space)

for (k in 1:n_evals) {
  start_time <- proc.time()  
  #simulate observations
  obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true,
                   sigma = sigma_true, Cx = Cx_true, Cy = Cy_true,
                   thetas = thetas_true, bd = bd)
  
  #find "best" center point 
  C_est <- best_C(C_poss, conts = obs$polys, thetas)
  
  #measure and store y
  #Make sets of lines, l, for  C_hat and get y's
  l_untrim <- make_l(C = C_est, thetas)
  l <- lapply(l_untrim, function(x){gIntersection(x, box)})
  l_lengths <- sapply(l, function(x){as.numeric(gLength(x))})
  
  #compute y's
  pts <- lapply(obs$polys, function(x) {pts_on_l(l, x, under = FALSE)})
  y <- sapply(pts, function(x){apply(x, 1, function(y){get_dist(y, C_est)})})
  
  #initial values for MCMC
  mu_ini <- apply(y, 1, mean)
  kappa_ini <- 1
  sigma_ini <-  apply(y, 1, sd)
  
  #non-scaler proposals for MCMC
  theta_dist_est <- theta_dist_mat(thetas)
  sigmaPropCov = sigmaProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)
  muPropCov <- muProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)
  
  #run MCMC
  fits <- ContouR::RunMCMC(nIter = n_iter, y,
                  mu = mu_ini, mu0 = mu0, Lambda0, muPropCov,
                  kappa = kappa_ini, betaKappa0, kappaPropSD,
                  sigma = sigma_ini, betaSigma0 = betaSigma0, sigmaPropCov,
                  gStart = g_start - 1, gEnd = g_end - 1,
                  thetaDist = theta_dist_est)

  #parameter estimates
  mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
  kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
  sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)

  #posterior field
  gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est,
                    sigma = sigma_est, Cx = C_est[1], Cy = C_est[2],
                    thetas = thetas, bd = bd)
  prob <- prob_field(gens$polys, nrows = n_grid, ncols = n_grid)
  
  #find credible intervals and compute coverage
  creds <- cred_regs(prob, cred_eval = cred_levels, nrows = n_grid, 
                     ncols = n_grid)
  if (k != 1) {
    cover <- cover + sapply(creds, 
                            function(x){eval_cred_reg(truth = test$polys[[k]],
                                                      cred_reg = x, 
                                                      center = c(Cx_true, Cy_true), 
                                                      p_test = p_test,
                                                      nrows = n_grid, ncols = n_grid)})
  } else {
    cover <- sapply(creds, 
                    function(x){eval_cred_reg(truth = test$polys[[k]],
                                              cred_reg = x, 
                                              center = c(Cx_true, Cy_true), 
                                              p_test = p_test,
                                              nrows = n_grid, ncols = n_grid,
                                              plotting = TRUE)})
  }
  end_time <- proc.time()
  elapse_time <- end_time  - start_time
  print(sprintf("Eval %i completed for task_id %i", k, task_id))
  print(elapse_time)
}

#save results
res_cover <- list("task" = task, "cover" = cover)
save(res_cover, 
     file = sprintf("/homes/direch/contours/Simulations/sim_results/nObsnGen/task%iShape%snObs%inGen%i.rda",
                    task_id, task$shape_name, n_obs, n_gen))
