#Test fits for a single set of observed contours
rm(list = ls())
set.seed(103)
library("sp")
library("fields")
library("rgeos")
library("Rcpp")
library("MASS")
source("R/gen_conts.R")
source("R/perpen_shift.R")
source("R/planes_intersections_polys.R")
source("R/misc.R")
source("R/credible_intervals.R")
sourceCpp("src/MCMC.cpp")


#truth
p_true <- 12
mu_true <- c(seq(.1, .4, length = p_true/4), seq(.4, .1, length = 3*p_true/4))
kappa_true <- 3
sigma_true <- c(seq(.005, .015, length = p_true/2),
                seq(.015, .005, length = p_true/2))
nu_true <- .01
Cx_true <- .5
Cy_true <- .5
theta_space <- 2*pi/p_true
theta1 <- theta_space/2
theta <- seq(theta1, 2*pi, by = theta_space)


#simulate observations
n_obs <- 25
obs <- gen_conts_simp(n_sim = n_obs, mu = mu_true, kappa = kappa_true, 
                     sigma = sigma_true, Cx = Cx_true, Cy = Cy_true,
                     theta1 = theta1)


#Priors (formal criteria needed)
mu0 <-  rep(.15, p_true); Lambda0 <- .05*diag(p_true) 
betaKappa0 <- 100
betaSigma0 <- .1
v0 <- .05
muC0 <- .5; tau20 <- .05^2
d0 <- .5

# #find observed intersection kernel
for (i in 1:n_obs) {
  kern_curr <- find_kernel(rbind(t(obs$coords[,,i]), t(obs$coords[,1,i])))
  if (i == 1) {
    kern <- kern_curr
  } else {
    kern <- gIntersection(kern, kern_curr)
  }
}
kern_pts <- t(kern@polygons[[1]]@Polygons[[1]]@coords)

C_ini <- gCentroid(kern)@coords


