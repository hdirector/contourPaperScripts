#Assess coverage under particular conditions
rm(list = ls())
set.seed(103)
library("ContouR")
start_time <- proc.time()

#Read in arguments
args <- commandArgs(trailingOnly =TRUE)
task_id <- as.numeric(args[1])
task_path <- args[2]
load(task_path)
task <- task_table[task_id,]

#simulation settings
n_obs <- task$n_obs
n_gen <- task$n_gen
cred_levels <- c(80, 90, 95)
n_cred <- length(cred_levels)
p_test <- task$p_test
n_evals <- task$n_evals
n_grid <- task$n_grid

#true parameters
pars <- get(as.character(task$shape_name))
p <- length(pars$mu)
mu_true <- pars$mu
p <- length(mu_true)
kappa_true <- pars$kappa
sigma_true <- pars$sigma
Cx_true <- pars$Cx
Cy_true <- pars$Cy
thetas_true <- pars$theta

#priors
mu0_indiv <- task$mu0
mu0 <- rep(mu0_indiv, p)
Lambda0_sigma2 <- task$Lambda0_sigma2
Lambda0 <- Lambda0_sigma2*diag(p) 
betaKappa0 <- task$betaKappa0
betaSigma0 <- task$betaSigma0

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
test <- gen_conts(n_sim = n_evals, mu = mu_true, kappa = kappa_true,
                  sigma = sigma_true, Cx = Cx_true, Cy = Cy_true,
                  thetas = thetas_true)

for (k in 1:n_evals) {  
  #simulate observations
  obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true,
                   sigma = sigma_true, Cx = Cx_true, Cy = Cy_true,
                   thetas = thetas_true)
  
  #find center point and associated angles
  C_est <- find_center(coords = obs$coords)
  thetas_all <- apply(obs$coords, 3, function(x){calc_angs(C = C_est,
                                                           coords = x)})
  stopifnot(max(apply(thetas_all, 1, sd)) <= 1e-10)
  thetas_est <- apply(thetas_all, 1, mean)
  
  #measure and store y
  temp <- XToWY(Cx = C_est[1], Cy = C_est[2], x = obs$coords, thetas_est)
  y <- temp$y
  
  #initial values for MCMC
  mu_ini <- apply(temp$y, 1, mean)
  kappa_ini = 5
  sigma_ini <-  apply(temp$y, 1, sd)
  
  #non-scaler proposals for MCMC
  theta_dist_est <- theta_dist_mat(thetas_est)
  sigmaPropCov = sigmaProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)
  muPropCov <- muProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)
  
  #run MCMC 
  fits <- RunMCMC(nIter = n_iter, y,
                  mu = mu_ini, mu0, Lambda0, muPropCov,
                  kappa = kappa_ini, betaKappa0, kappaPropSD,
                  sigma = sigma_ini, betaSigma0, sigmaPropCov,
                  gStart = g_start - 1, gEnd = g_end - 1, 
                  thetaDist = theta_dist_est)
  
  #parameter estimates
  mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
  kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
  sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)

  #posterior field
  gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est,
                    sigma = sigma_est, Cx = C_est[1], Cy = C_est[2],
                    thetas = thetas_est)
  prob <- prob_field(gens$polys, nrows = n_grid, ncols = n_grid)
  
  #find credible intervals and compute coverage
  creds <- cred_regs(prob, cred_levels, nrows = n_grid, ncols = n_grid)
  if (k != 1) {
    cover <- cover + sapply(creds, 
                            function(x){eval_cred_reg(truth = test$polys[[k]],
                                                      cred_reg = x, 
                                                      center = c(Cx_true, Cy_true), 
                                                      p_test = p_test)})
  } else {
    cover <- sapply(creds, 
                    function(x){eval_cred_reg(truth = test$polys[[k]],
                                              cred_reg = x, 
                                              center = c(Cx_true, Cy_true), 
                                              p_test = p_test)})
  }
  print(sprintf("Eval %i completed for task_id %i", k, task_id))
}
end_time <- proc.time()
end_time  - start_time
print(sprintf("task_id %i completed", task_id))
print(elapse_time)
#save(cover, file = "")
