#Test Script shows that MCMC works with large number of samples and fixed C 
#and theta

#Assess coverage under particular conditions
rm(list = ls())
set.seed(103)
library("ContouR")
library("coda")

#Read in simulation settings
task_id <- 4
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/nObsnGenTaskTable.rda"
load(task_path)
task <- task_table[task_id,]
print(sprintf("task_id: %i", task_id))
attach(task)
n_iter <- 50000
burn_in <- 25000
n_obs <- 100

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

test <- gen_conts(n_sim = n_evals, mu = mu_true, kappa = kappa_true,
                  sigma = sigma_true, C = C_true, thetas = thetas_true, 
                  bd = bd)

#use the correct theta and p to confirm that sampler works 
theta_space_true <- 2*pi/p_true
thetas_true <- seq(theta_space_true/2, 2*pi, theta_space_true)

start_time <- proc.time()  
#simulate observations
obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true,
                 sigma = sigma_true, C = C_true, thetas = thetas_true,
                 bd = bd)

#measure and store y
#Make sets of lines, l, for  C_hat and get y's
l_untrim <- make_l(C = C_true, thetas_true)
l <- lapply(l_untrim, function(x){gIntersection(x, box)})
l_lengths <- sapply(l, function(x){as.numeric(gLength(x))})

#compute y's
pts <- lapply(obs$polys, function(x) {pts_on_l(l, x, under = FALSE)})
y <- sapply(pts, function(x){apply(x, 1, function(y){get_dist(y, C_true)})})

#initial values for MCMC
mu_ini <- apply(y, 1, mean)
kappa_ini <- 1
sigma_ini <-  apply(y, 1, sd)

#priors
mu0_indiv <- task$mu0
mu0 <- rep(mu0_indiv, p_true)
Lambda0 <- Lambda0_sigma2*diag(p_true) 
betaSigma0 <- rep(task$betaSigma0, p_true)


#non-scaler proposals for MCMC
theta_dist_est <- theta_dist_mat(thetas_true)
sigmaPropSD = sqrt(diag(sigmaProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)))
muPropSD <- sqrt(diag(muProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)))



#run MCMC
start_time <- proc.time()
fits <- RunMCMC(nIter = n_iter, y,
                mu = mu_ini, mu0 = mu0, Lambda0, muPropSD,
                kappa = kappa_ini, betaKappa0, kappaPropSD = .05,
                sigma = sigma_ini, betaSigma0 = betaSigma0, sigmaPropSD,
                thetaDist = theta_dist_est)
end_time <- proc.time()
run_time <- end_time - start_time
print(sprintf("fitted eval %i, time: %s", k, run_time[3]))


#parameter estimates
mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)

par(mfrow = c(1, 2))
Sigma_obs <- cov(t(y))
Sigma_ini <- compSigma(sigma_ini, kappa_ini, theta_dist_est)
Sigma_est <- compSigma(sigma_est, kappa_est, theta_dist_est)
image.plot(Sigma_est, zlim  = c(0, .008))
image.plot(Sigma_obs, zlim = c(0, .008))

par(mfrow = c(1, 1))
plot(sigma_est, type= "l")
points(sigma_ini, type= "l", col = "blue")
points(sigma_true, type= "l", col = 'red')
points(betaSigma0, type = "l", lty = 2)
plot(mu_est, type= "l")
points(mu_ini, col= 'blue', type= "l")
plot(fits$kappa, type= "l")
kappa_est; kappa_true


