#Simulation study, increasing number of contours
rm(list = ls())
set.seed(103)

comp_coverage <- function(settings, pars, p_test, n_evals) {
  library("PolygonKernels")
  #varied settings
  n_obs <- settings[1]
  n_gen <- settings[2]

  #true parameters
  mu_true <- pars$mu
  p_true <- length(mu_true)
  kappa_true <- pars$kappa
  sigma_true <- pars$sigma
  Cx_true <- pars$Cx
  Cy_true <- pars$Cy
  
  theta_space <- 2*pi/p_true
  theta1_true <- theta_space/2
  thetas_true <- seq(theta1_true, 2*pi, by = theta_space)
  
  #Priors (formal criteria needed)
  p_est <- p_true
  mu0 <-  rep(.15, p_est); Lambda0 <- .05*diag(p_est) 
  betaKappa0 <- 100
  betaSigma0 <- .1

  #Fixed simulation info 
  n_iter <- 30000
  burn_in <- 20000
  g_space <- 1
  g_start <- seq(1, p_est, by = g_space)
  g_end <- c(seq(g_space, p_true, by = g_space), p_est)
  n_cred <- length(cred_eval)
  n_grid <- 300
  
  #sample "true" contour
  test <- gen_conts(n_sim = n_evals, mu = mu_true, kappa = kappa_true,
                    sigma = sigma_true, Cx = Cx_true, Cy = Cy_true,
                    theta1 = theta1_true)
  
  
  for (k in 1:n_evals) {  
    #simulate observations
    obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true,
                     sigma = sigma_true, Cx = Cx_true, Cy = Cy_true,
                     theta1 = theta1_true)

    #find center point and associated angles
    C_est <- find_center(obs_coords = obs$coords)
    print(C_est)
    thetas_all <- apply(obs$coords, 3, function(x){calc_angs(C = C_est,
                                                             coords = x)})
    stopifnot(max(apply(thetas_all, 1, sd)) <= .1)
    thetas_est <- apply(thetas_all, 1, mean)
    x = array(sapply(obs$polys, function(x){paral_pts(p = p_est, poly = x, 
                                                      center = C_est,  r = 10)}),
              dim = c(2, p_est, n_obs))
    
    #initial values
    temp <- XToWY(Cx = C_est[1], Cy = C_est[2], x = x, thetas_est)
    mu_ini <- rep(mean(apply(temp$y, 1, mean)), p_est) 
    kappa_ini = 5
    sigma_ini <-  rep(mean(apply(temp$y, 1, sd)), p_est)

    #non-scaler proposals
    theta_dist_est <- compThetaDist(p_est, theta_space)
    sigmaPropCov = .005*compSigma(sigma_ini, kappa_ini, theta_dist_est)
    muPropCov <- .001*compSigma(sigma_ini, kappa_ini, theta_dist_est)
    
    ####TO DO: CHANGE COMP THETA DIST IN MCMC CODE TO ALLOW FOR VARIED ANGLE SPACING
    fits <- RunMCMC(nIter = n_iter, x = x,
                    mu = mu_ini, mu0, Lambda0, muPropCov,
                    kappa = kappa_ini, betaKappa0, kappaPropSD = .05,
                    sigma = sigma_ini, betaSigma0, sigmaPropCov,
                    Cx = C_est[1], Cy = C_est[2],
                    theta1 = thetas_est[1],
                    gStart = g_start - 1, gEnd = g_end - 1)
    
    #parameter estimates
    mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
    kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
    sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)
    rm(fits) #save memory
    
    #posterior field
    gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est,
                      sigma = sigma_est, Cx = Cx_true, Cy = Cy_true,
                      theta1 = theta1_true)
    prob <- prob_field(gens$polys, nrows = n_grid, ncols = n_grid)
    rm(gens) #save memory
    

    #result matrix
    creds <- cred_regs(prob, cred_eval, nrows = n_grid, ncols = n_grid)
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
    print(sprintf("n_eval: %i, n_obs: %i, n_gen %i", k, n_obs, n_gen))
  }
  return(cover)
}

#set up varied simulation parameters to test
n_obs_poss <- c(5, 20, 40)
n_gen_poss <- 100#c(20, 50, 100)
settings_list <- list()
settings_tab <- expand.grid(n_obs_poss, n_gen_poss)
for (i in 1:nrow(settings_tab)) {
  settings_list[[i]] <- as.numeric(settings_tab[i,])
}

#run cases in parallel
cred_eval <- c(80, 90, 95)
p_test <- 24
n_evals <- 100


######single test case######
# test <- comp_coverage(settings = settings_list[[1]], pars = stop_sign, p_test = p_test, 
#                        n_evals = n_evals)
#########################

##########################
#Real evaluation 
###########################
library("parallel")
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores, outfile = "")
clusterExport(cl, c("comp_coverage", "p_test", "n_evals", "cred_eval"))
res_stop <- parLapply(cl, settings_list, function(x){comp_coverage(settings = x, 
                                                         pars = stop_sign,
                                                         p_test = p_test,
                                                         n_evals = n_evals)})
save(res_stop, file  = "~/desktop/res_stop_sign.rda")
res_tie <- parLapply(cl, settings_list, function(x){comp_coverage(settings = x, 
                                                              pars = tie,
                                                              p_test = p_test,
                                                              n_evals = n_evals)})
save(res_tie, file  = "~/desktop/res_tie.rda")
stopCluster(cl)

#evaluate results
res_sum <- t(sapply(res, function(x){apply(x, 2, function(y){mean(y/n_evals)})}))
library("tidyverse")
res_df <- data.frame("n_obs" = rep(as.factor(par_tab[,1]), 3),
                      "n_gen" = rep(as.factor(par_tab[,2]), 3),
                      "nom_cover" = rep(cred_eval, each = nrow(par_tab)),
                      "cover" = 100*as.vector(res_sum[,1:3]),
                      "cat" = rep("data", nrow(par_tab)))

temp <- filter(res_df, n_gen == 20)
xtable(temp[,1:4], digits = 1)

#pdf("/users/hdirector/desktop/obsExper.pdf")
ggplot(filter(res_df, n_gen == 100),
       aes(x = nom_cover, y = jitter(cover), group = n_obs, col = n_obs)) +
  geom_point() +
  ggtitle("Coverage Improves With Sample Size (Jittered")
#dev.off()

#pdf("/users/hdirector/desktop/simExper.pdf")
ggplot(filter(res_df, n_obs == 50),
       aes(x = nom_cover, y = jitter(cover), group = n_gen, col = n_gen)) +
  geom_point() +
  ggtitle("Coverage Improves With More Simulations")
#dev.off()


