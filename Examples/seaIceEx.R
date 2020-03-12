#Assess coverage for sea ice example
rm(list = ls())
set.seed(103)
library("IceCast")
library("ContouR")
library("raster")
library("coda")
library("fields")
library("viridis")

#data about original grid
xmn = -3850; xmx = 3750; ymn = -5350; ymx = 5850
box_size <- 25

eps <- .1 #distance away from edge of box for rescale
start_year <- 2010
end_year <- 2017
n_grid <- 300
p_init <- 200
p_eval <- 100
err_prop <- .05

#derived settings
years <- start_year:end_year
n_years <- length(years)
n_train <- n_years - 1
box <- bbox(eps) 

#make a grid of points at the center of each square 
x_cent <- seq(xmn + 13.5, xmx, 150)
y_cent <- seq(ymn + 13.5, ymx, 150)
grid <- expand.grid(x_cent, y_cent)

#regions in and around central Arctic
reg <- reg_info$regions[[1]]
reg@polygons[[1]]@ID <- "CentArc"
reg <- spRbind(reg, reg_info$regions[[2]])
reg <- rm_holes(raster::aggregate(reg))
reg_coords <- reg@polygons[[1]]@Polygons[[1]]@coords

#load data
all <- IceCast::read_monthly_BS(start_year, end_year, 
                                "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/bootstrapV3_1/",
                                version = 3.1)
sep <- all[,9,,]

coords <- list()
for (i in 1:n_years) {
  #find coordinates of largest contigous region in Central Arctic
  temp <- disaggregate(get_region(sep[i,,], dat_type = "bootstrap", level = 15))
  temp <- temp[which.max(gArea(temp, byid = TRUE))]
  temp <- disaggregate(rm_holes(keep_poly(gIntersection(reg, temp))))
  temp <- temp[which.max(gArea(temp, byid = T))]
  coords[[i]] <- temp@polygons[[1]]@Polygons[[1]]@coords
}

#rescale everything 
temp_rescale <- rescale(coords = coords, eps = eps, box_size = box_size, 
                        grid = grid, bd = reg_coords)

coords_scale <- temp_rescale$coords_scale
grid_scale <- temp_rescale$grid_scale
bd_scale <- temp_rescale$bd_scale
bd_scale <- bd_scale[nrow(bd_scale):1, ] #reverse order to match generated contours
box_size <- temp_rescale$box_size

#Switch to Spatial objects
grid_pts_all <- SpatialPoints(grid_scale)
conts <- lapply(coords_scale, function(x){make_poly(x, "contours")})

#identify section with no variability and remove from 
bd_poly <- make_poly(bd_scale, "bd_poly")
no_var <- gIntersection(as(bd_poly, "SpatialLines"), 
                        as(conts[[1]], "SpatialLines"))
for (i in 2:n_years) {
  no_var <- gIntersection(no_var, conts[[i]])
}

#Find C and thetas for evaluation contour
theta_space_eval <- 2*pi/p_eval
thetas_eval_all <- seq(theta_space_eval/2, 2*pi, by = theta_space_eval)
area_tol <- err_prop*mean(sapply(conts, gArea)) 
C_eval <- best_C(bd = bd_scale, conts = conts, thetas = thetas_eval_all,
                 area_tol = area_tol)
l_untrim_eval <- make_l(C = C_eval, thetas_eval_all)
l_eval <- lapply(l_untrim_eval, function(x){gIntersection(x, bd_poly)})
on_can_arch_eval <- sapply(l_eval, function(x){gIntersects(no_var, x)})
l_eval <- l_eval[!on_can_arch_eval]
thetas_eval <- thetas_eval_all[!on_can_arch_eval]

#input settings
for (test_year in 2010:2017)  {
  eval_ind <- which(years == test_year)
  theta_space_init <- 2*pi/p_init
  thetas_init <- seq(theta_space_init/2, 2*pi, by = theta_space_init)
  
  #Seperate out training and testing set
  train <- conts
  train[[eval_ind]] <- NULL
  eval <- conts[[eval_ind]]
  
  #find estimated center point, p, and theta
  area_tol <- err_prop*mean(sapply(train, gArea)) 
  C_est <- best_C(bd = bd_scale, conts = train, thetas = thetas_init, 
                  area_tol = area_tol)
  p_est_all <- reduce_p(C = C_est, conts = train, area_tol = area_tol, 
                    p = p_init, red_prop = .05)
  theta_space_est <- 2*pi/p_est_all
  thetas_est_all <- seq(theta_space_est/2, 2*pi, theta_space_est)

  #Make sets of lines, l, for  C_est and get y's
  l_untrim <- make_l(C = C_est, thetas_est_all)
  l <- lapply(l_untrim, function(x){gIntersection(x, bd_poly)})
  on_can_arch <- sapply(l, function(x){gIntersects(no_var, x)})
  l_lengths_all <- sapply(l, function(x){as.numeric(gLength(x))})
  l_rm <- which(on_can_arch)
  n_rm <- length(l_rm)
  l_keep <- which(!on_can_arch)
  l <- l[l_keep]
  
  #adjust parameters for only those that will be fit
  l_lengths <- l_lengths_all[l_keep]
  thetas_est <- thetas_est_all[l_keep]
  p_est <- p_est_all - n_rm
  
  #priors
  mu0 <- (2/3)*l_lengths
  Lambda0_sigma2 <- .5^2
  Lambda0 <- Lambda0_sigma2*diag(p_est) 
  betaKappa0 <- 10
  betaSigma0 <- (l_lengths/2)/qnorm(.995)
  
  #compute y's
  pts <- lapply(train, function(x) {pts_on_l(l, x, under = FALSE)})
  y <- sapply(pts, function(x){apply(x, 1, function(y){get_dist(y, C_est)})})
  
  #initial values for MCMC
  mu_ini <- apply(y, 1, mean)
  kappa_ini <- .5
  sigma_obs <-  apply(y, 1, sd)
  sigma_ini <- sigma_obs
  sigma_ini[sigma_ini == 0] <- .0001

  #MCMC settings
  n_iter <- 50000
  burn_in <- 25000
  sigmaProp_sigma2 <- .05
  muProp_sigma2 <- .01
  kappaPropSD <- .01
  cz <- contig_zero(sigma_obs)
  g_space <- 1
  g_start <- unlist(apply(cz, 1, function(x){seq(x[1], x[2], g_space)}))
  g_end <- c(g_start[2:length(g_start)] - 1, p_est)
  
  #non-scaler proposals for MCMC
  theta_dist_est <- theta_dist_mat(thetas_est)
  sigmaPropCov = sigmaProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)
  muPropCov <- muProp_sigma2*compSigma(sigma_ini, kappa_ini, theta_dist_est)
  
  #run MCMC
  start_time <- proc.time()
  fits <- ContouR::RunMCMC(nIter = n_iter, y,
                          mu = mu_ini, mu0 = mu0, Lambda0, muPropCov,
                          kappa = kappa_ini, betaKappa0, kappaPropSD,
                          sigma = sigma_ini, betaSigma0 = betaSigma0, sigmaPropCov,
                          gStart = g_start - 1, gEnd = g_end - 1,
                          thetaDist = theta_dist_est)
  end_time <- proc.time()
  run_time <- end_time - start_time
  print(sprintf("fitted year %i complete, time: %s", test_year, run_time[3]))
  
  
  #parameter estimates
  mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
  kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
  sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)
  
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
  # 
  
  #########################
  #Evaluate fits
  #########################
  Sigma_obs <- cov(t(y))
  Sigma_ini <- compSigma(sigma_ini, kappa_ini, theta_dist_est)
  Sigma_est <- compSigma(sigma_est, kappa_est, theta_dist_est)
  # par(mfrow = c(1, 2))
  # image.plot(Sigma_est)
  # image.plot(Sigma_obs)
# 
#   plot(sigma_est, type= "l")
#   points(sigma_ini, type= "l", col = "blue")
#   points(betaSigma0, type = "l", lty = 2)
#   plot(mu_est, type= "l")
#   points(mu_ini, col= 'blue', type= "l")

  #put back lines along Canadian archipelago
  mu_est_all <- c(mu_est[1:(l_rm[1] - 1)], l_lengths_all[l_rm],
                  mu_est[l_rm[1]:p_est])
  sigma_est_all <- c(sigma_est[1:(l_rm[1] - 1)], rep(0, n_rm),
                 sigma_est[l_rm[1]:p_est])
  
  #posterior field
  n_gen <- 100
  gens <- gen_conts(n_sim = n_gen, mu = mu_est_all, kappa = kappa_est, 
                    sigma = sigma_est_all, C = C_est, thetas = thetas_est_all,
                    bd = bd_scale)
  

  #set up grid to align with data
  grid_need <- 1/box_size
  xmn <- ymn <- 0 
  xmx <- ymx <- 1 
  x_bds <- seq(xmn, xmx, by = box_size[1])
  y_bds <- seq(ymn, ymx, by = box_size[2])

  prob <- prob_field(polys = gens$polys, nrows = length(y_bds) - 1, 
                     ncols = length(x_bds) - 1,  xmn = xmn, xmx = xmx, 
                     ymn = ymn, ymx = ymx)

  #compare results to generated contours
  pdf(sprintf("/users/hdirector/desktop/probAndGen%i.pdf", test_year),
      width = 8.5, height = 4)
  par(mfrow = c(1, 2))
  image.plot(x_bds, y_bds, prob, xaxt = "n", yaxt = "n", col = viridis(20),
             main = "Probability Field")
  for (i in 1:n_train) {
    plot(train[[i]], add = T, border = "black")
  }
  plot(eval, add = T, border = "red")
  
  
  creds <- cred_regs(prob, cred_eval = c(80, 90, 95), nrows = n_grid, ncols = n_grid)
  plot(creds[[3]], col = 'navy', lwd = 2, border = "navy", 
       main = "Credible Intervals")
  plot(creds[[2]], col = 'blue', lwd = 2, add = TRUE, border = "blue")
  plot(creds[[1]], col = 'lightblue', lwd = 2, add = TRUE, border = "lightblue")
  
  for (i in 1:n_train) {
    plot(train[[i]], add = TRUE, border = "black")
  }
  plot(eval, add = TRUE, border = 'red')
  legend("topright",  fill = c("navy", "blue", "lightblue", "black", "red"),
         legend = c("80%","90%", "95", "Training Contours", "Test Contour"), 
         cex = .45)
  dev.off()
  
  #####################
  #formal evaluation
  #####################
  #evaluate coverage
  pdf(file = sprintf("/users/hdirector/desktop/cred_ints%i.pdf", test_year),
              width = 8.5, height = 4)
  par(mfrow = c(1, 3))
  cover <- sapply(creds, 
         function(x){eval_cred_reg(truth = eval, cred_reg = x, 
                                   center = C_eval,  thetas = thetas_eval, 
                                   nrows = n_grid, ncols = n_grid,
                                   plotting = TRUE)})
  mtext(text = "80% Credible Interval", side = 3, line = -3, at = .18,
        outer = TRUE)
  mtext(text = "90% Credible Interval", side = 3, line = -3, at = .5, 
        outer = TRUE)
  mtext(text = "95% Credible Interval", side = 3, line = -3, at = .82,
        outer = TRUE)
  dev.off()
  
  save(cover, file = sprintf("/users/hdirector/desktop/cover_year%i.rda", 
                             test_year))
  print(sprintf("completed analysis of year %i", test_year))
}

# pdf("/users/hdirector/Dropbox/Contours/contourPaperScripts/figures/sea_ice_field.pdf",
#     height = 4, width = 4)
# par(mfrow = c(1, 1), mar = rep(1, 4), oma = c(0, 0, 0, 1))
# image.plot(x_bds, y_bds, prob, xaxt = "n", yaxt = "n", col = viridis(20),
#            xlim = c(0.1, 0.9), ylim = c(0.1, 0.9), xaxt = "n", yaxt = "n",
#            xlab = "", ylab = "")
# plot(eval, add = T, border = "red")
# outside <- gDifference(box, bd_poly)
# plot(outside, add = T, col = 'beige')
# plot(box, add = T, border = 'black')
# dev.off()
# 
