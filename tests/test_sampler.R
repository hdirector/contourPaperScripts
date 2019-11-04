#Test fits for a single set of observed contours
rm(list = ls())
set.seed(103)
library("PolygonKernels")

#truth
p_true <- 24
mu_true <- c(seq(.1, .4, length = p_true/4), seq(.4, .1, length = 3*p_true/4))
kappa_true <- 3
sigma_true <- c(seq(.005, .015, length = p_true/2),
                seq(.015, .005, length = p_true/2))
Cx_true <- .5
Cy_true <- .5
theta_space <- 2*pi/p_true
theta1_true <- theta_space/2
thetas_true <- seq(theta1_true, 2*pi, by = theta_space)

#simulate observations
n_obs <- 25
obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true, 
                 sigma = sigma_true, Cx = Cx_true, Cy = Cy_true,
                 theta1 = theta1_true)

#pdf("Figures/test_sampler/test_sampler_obs.pdf")
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", col = 'white',
     main = "Observed Contours")
plot(obs$polys[[1]], add = T, lwd = .25)
for (i in 2:n_obs) {
  plot(obs$polys[[i]], add = T, lwd = .25)
}
points(Cx_true, Cy_true, col = 'blue')
points(C_est[1], C_est[2], col = 'red')
#dev.off()

#Priors (formal criteria needed)
p_est <- 24
mu0 <-  rep(.15, p_est); Lambda0 <- .05*diag(p_est) 
betaKappa0 <- 100
betaSigma0 <- .1
muC0 <- .2; sigmaC0 <- .01

#find center point and associated angles
C_est <- find_center(obs$coords)
thetas_all <- apply(obs$coords, 3, function(x){calc_angs(C = C_est, coords = x)})
stopifnot(max(apply(thetas_all, 1, sd)) <= .001)
thetas_est <- thetas_all[,1]
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

#Fixed simulation info 
n_iter <- 30000
burn_in <- 20000
g_space <- 1
g_start <- seq(1, p_est, by = g_space)
g_end <- c(seq(g_space, p_true, by = g_space), p_est)


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


#mu evaluation
for (i in 1:p_est) {
  plot(fits$mu[i,], type = "l")
  abline(h = mu_true[i], col = 'red')
}
#png(file = "Figures/test_sampler/test_sampler_mu2_tp.png")
plot(fits$mu[2,], type = "l", xlab = "Iteration",  ylab = expression(mu[2])) 
abline(h = mu_true[2], col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
#dev.off()
plot(thetas_true, mu_true, col = 'blue')
points(thetas_est, mu_est, pch = 20)
fits$muRate

#sigma evaluation
for (i in 1:p_true) {
  plot(fits$sigma[i,], type= "l")
  abline(h = sigma_true[i], col = 'red')
}
#png(file = "Figures/test_sampler/test_sampler_sigma2_tp.png")
plot(fits$sigma[2,], type = "l", xlab = "Iteration",  ylab = expression(sigma[2])) 
abline(h = sigma_true[2], col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
#dev.off()
plot(thetas_true, sigma_true, col = 'blue', ylim = c(.004, .018))
points(thetas_est, sigma_est, pch = 20)
fits$sigmaRate

#kappa evaluation
#png(file = "Figures/test_sampler/test_sampler_kappa_tp.png")
plot(fits$kappa, type = "l", xlab = "Iteration",  ylab = expression(kappa)) 
abline(h = kappa_true, col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
#dev.off()
plot(fits$kappa, type = "l")
abline(h = kappa_true, col = 'red')
fits$kappaRate


#posterior field
n_gen <- 200
gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est,
                  sigma = sigma_est, Cx = C_est[1], Cy = C_est[2],
                  theta1 = thetas_est[1])
prob <- prob_field(polys = gens$polys, nrows = 300, ncols = 300)


#pdf("figures/test_sampler/test_sampler_res.pdf", height = 4, width = 8)
par(mfrow = c(1, 2))
image.plot(prob)
#plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", col = 'white')
plot(gens$polys[[1]], add = T, lwd = .5, border = 'green')
plot(obs$polys[[1]], add = T, lwd = .5)
for (i in 2:25) {
  #plot(gens$polys[[i]], add = T, lwd = .5, border = 'green')
  plot(obs$polys[[i]], add = T, lwd = .5)
}
legend("bottom", horiz = T, fill = c("green", "black"),
       legend = c("'Observed'", "Generated"), cex = .75)
#dev.off()

