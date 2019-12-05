#Assess coverage under particular conditions
rm(list = ls())
set.seed(103)
library("ContouR")
library("IceCast")
library("raster")

#settings
p <- 100 #number of lines on which to build
eps <- .1 #distance away from edge of box for rescale
box <- bbox(eps) 
theta_dist <- 2*pi/p
thetas <- seq(theta_dist/2, 2*pi, by = theta_dist)

#1:10 = 1980:1989, 11:20 = 1990:1999, 21:30 = 2000:2009, 31:38 = 2010:2017 
year_ind <- 2010:2017


#load data
all <- IceCast::read_monthly_BS(2010, 2017, 
                                "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/bootstrapV3_1/",
                                version = 3.1)
sep <- all[,9,,]
n_obs <- dim(sep)[1]

#regions in and around central Arctic
reg <- reg_info$regions[[1]]
reg@polygons[[1]]@ID <- "CentArc"
reg <- spRbind(reg, reg_info$regions[[2]])
reg <- rm_holes(raster::aggregate(reg))

#make a grid of points at the center of each square (
xmn = -3850; xmx = 3750; ymn = -5350; ymx = 5850
x_cent <- seq(xmn + 13.5, xmx, 150)
y_cent <- seq(ymn + 13.5, ymx, 150)
grid <- expand.grid(x_cent, y_cent)

coords <- list()
for (i in 1:n_obs) {
  #find coordinates of largest contigous region in Central Arctic
  temp <- disaggregate(get_region(sep[i,,], dat_type = "bootstrap", level = 15))
  temp <- temp[which.max(gArea(temp, byid = TRUE))]
  temp <- disaggregate(rm_holes(keep_poly(gIntersection(reg, temp))))
  temp <- temp[which.max(gArea(temp, byid = T))]
  coords[[i]] <- temp@polygons[[1]]@Polygons[[1]]@coords
}

#rescale everything 
temp_rescale <- rescale(coords, eps, grid)
coords_scale <- temp_rescale$coords_scale
grid_scale <- temp_rescale$grid_scale

#Switch to Spatial objects
grid_pts_temp <- SpatialPoints(grid_scale)
conts <- lapply(coords_scale, function(x){make_poly(x, "contours")})

#identify set of test points that are in every contour
test <- sapply(conts, function(x){gIntersects(grid_pts_temp, x, byid = TRUE)})
keep <- apply(test, 1, function(x){all(x)})
grid_pts <- SpatialPoints(grid_scale[keep,])
n_test_pts <- nrow(grid_pts@coords)

#Compute areas in error (area_out)
#TO DO: update to use function y_on_l
# C_test <- grid_pts@coords
# area_out <- matrix(nrow = n_test_pts, ncol = n_obs)
# for (i in 1:n_test_pts) {
#    l <- make_l(C = C_test, theta = theta)
#   for (j in 1:n_obs) {
#     on_l <- lapply(l, function(x){gIntersection(conts[[j]], x)})
#     pts_on_l <- t(sapply(on_l, function(x){x@lines[[1]]@Lines[[1]]@coords[2,]}))
#     test_cont <- make_poly(pts_on_l, "test_cont")
#     diff_reg1 <- gDifference(test_cont, conts[[j]])
#     diff_reg2 <- gDifference(conts[[j]], test_cont)
#     area_out[i, j] <- gArea(diff_reg1) + gArea(diff_reg2)
#     print(c(i, j))
#   }
# }
load("/users/hdirector/desktop/area_out.rda")

#determine center ponts for each decade
max_area_out <- apply(area_out[,31:38], 1, max) #Re-run area_out
opt_ind <- which.min(max_area_out)
opt_areas <- area_out[opt_ind, 31:38]
C_hat <- grid_pts@coords[opt_ind,]

#Make sets of lines, l, for each C_hat and get y's
l_untrim <- make_l(C = C_hat, thetas)
l <- lapply(l_untrim, function(x){gIntersection(x, box)})
l_lengths <- sapply(l, function(x){as.numeric(gLength(x))})

# #Figure showing that data are approximately star-shaped
# pdf("/users/hdirector/desktop/approxPerform.pdf")
# par(mfrow = c(2, 3), oma = c(0, 1, 0, 0), mar = c(0, 4, 1, 0))
# for (i in 1:n_case) {
#    l <- make_l(C = C_hat[i,], theta)
#   #plot 3 cases (min, max, and, one in the middle)
#   n_inds <- length(inds_i)
#   for (j in c(1, floor(n_inds/2), n_inds)) {
#     curr <- which(order(opt_areas_i) == j)
#     cont_ind <- dec_ind[i, 1] + curr - 1
#     on_l <- lapply(l, function(x){gIntersection(conts[[cont_ind]], x)})
#     pts_on_l <- t(sapply(on_l, function(x){x@lines[[1]]@Lines[[1]]@coords[2,]}))
#     test_cont <- make_poly(pts_on_l, "test_cont")
#     if (j == 1) {
#       plot(conts[[cont_ind]], main = sprintf("%i-th Ranked Area", j), 
#            xlim = c(.1, .9), ylim = c(.1, .9), lwd = 1, ylab = dec_names[i])
#     } else {
#       plot(conts[[cont_ind]], main = sprintf("%i-th Ranked Area", j), 
#            xlim = c(.1, .9), ylim = c(.1, .9), lwd = 1)
#     }
#     plot(test_cont, add = T, border = 'blue', lwd = 1)
#     points(pts_on_l, col = 'blue', pch = 20)
#     points(C_hat[i, 1], C_hat[i, 2], col = 'red',pch = 20)
#   }
# }
# dev.off()

#priors
mu0 <- l_lengths/2
Lambda0_sigma2 <- .05
Lambda0 <- Lambda0_sigma2*diag(p) 
betaKappa0 <- 10
betaSigma0 <- (l_lengths/2)/qnorm(.995)

#compute y's
pts <- lapply(conts, function(x) {pts_on_l(l, x)})
y <- sapply(pts, function(x){apply(x, 1, function(y){get_dist(y, C_hat)})})


#initial values for MCMC
mu_ini <- apply(y, 1, mean)
kappa_ini <- 2.5
sigma_obs <-  apply(y, 1, sd)
sigma_ini <- sigma_obs
sigma_ini[sigma_ini == 0] <- .0001

#MCMC settings
n_iter <- 10000
burn_in <- 5000
sigmaProp_sigma2 <- .1
muProp_sigma2 <- .01
kappaPropSD <- .01
cz <- contig_zero(sigma_obs)
g_space <- 5
g_start <- unlist(apply(cz, 1, function(x){seq(x[1], x[2], g_space)}))
g_end <- c(g_start[2:length(g_start)] - 1, p)

#non-scaler proposals for MCMC
theta_dist_est <- theta_dist_mat(thetas)
sigmaPropCov = sigmaProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)
muPropCov <- muProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)


#run MCMC
start_time <- proc.time()
fits <- RunMCMC(nIter = n_iter, y,
                mu = mu_ini, mu0 = mu0, Lambda0, muPropCov,
                kappa = kappa_ini, betaKappa0, kappaPropSD,
                sigma = sigma_ini, betaSigma0 = betaSigma0, sigmaPropCov,
                gStart = g_start - 1, gEnd = g_end - 1,
                thetaDist = theta_dist_est, alpha = 0.5)
end_time <- proc.time()

#parameter estimates
mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)

Sigma_obs <- cov(t(y))
Sigma_ini <- compSigma(sigma_ini, kappa_ini, alpha = .5, theta_dist_est)
Sigma_est <- compSigma(sigma_est, kappa_est, alpha = .5, theta_dist_est)
par(mfrow = c(1, 2))
image.plot(Sigma_est)
image.plot(Sigma_obs)


plot(sigma_est, type= "l")
points(sigma_ini, type= "l", col = "blue")
points(betaSigma0, type = "l", lty = 2)
plot(mu_est, type= "l")
points(mu_ini, col= 'blue', type= "l") 

plot(fits$kappa, type= "l")
plot(fits$sigma[100,], type= "l")

#posterior field
n_gen <- 100
gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est,
                  sigma = sigma_est, Cx = opt_C[1], Cy = opt_C[2],
                  thetas = thetas)
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


#save results
res_cover <- list("task" = task, "cover" = cover)
save(res_cover, 
     file = sprintf("/homes/direch/contours/Simulations/sim_results/nObsnGen/task%iShape%snObs%inGen%i.rda",
                    task_id, task$shape_name, n_obs, n_gen))
