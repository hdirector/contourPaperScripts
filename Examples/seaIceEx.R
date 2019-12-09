#Assess coverage under particular conditions
rm(list = ls())
set.seed(103)
library("ContouR")
library("IceCast")
library("raster")
library("coda")
library("fields")
library("viridis")

#input settings
for (test_year in 2010:2017)  {
  p <- 100 #number of lines on which to build
  eps <- .1 #distance away from edge of box for rescale
  start_year <- 2010
  end_year <- 2017
  n_grid <- 300
  p_test <- 50
  
  #derived settings
  years <- start_year:end_year
  n_years <- length(years)
  eval_ind <- which(years == test_year)
  train_ind <- 1:n_years
  train_ind <- train_ind[-eval_ind]
  box <- bbox(eps) 
  theta_dist <- 2*pi/p
  thetas <- seq(theta_dist/2, 2*pi, by = theta_dist)
  
  #load data
  all <- IceCast::read_monthly_BS(start_year, end_year, 
                                  "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/bootstrapV3_1/",
                                  version = 3.1)
  sep <- all[,9,,]
  n_train <- n_years - 1
  
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
  for (i in 1:n_years) {
    #find coordinates of largest contigous region in Central Arctic
    temp <- disaggregate(get_region(sep[i,,], dat_type = "bootstrap", level = 15))
    temp <- temp[which.max(gArea(temp, byid = TRUE))]
    temp <- disaggregate(rm_holes(keep_poly(gIntersection(reg, temp))))
    temp <- temp[which.max(gArea(temp, byid = T))]
    coords[[i]] <- temp@polygons[[1]]@Polygons[[1]]@coords
  }
  
  #rescale everything 
  temp_rescale <- rescale(coords = coords, eps = eps, grid = grid)
  coords_scale <- temp_rescale$coords_scale
  grid_scale <- temp_rescale$grid_scale
  
  #Switch to Spatial objects
  grid_pts_all <- SpatialPoints(grid_scale)
  conts <- lapply(coords_scale, function(x){make_poly(x, "contours")})
  
  #Seperate out training and testing set
  train <- conts
  train[[eval_ind]] <- NULL
  eval <- conts[[eval_ind]]
  
  #identify set of test points that are in every contour
  poss_C_pts <- sapply(train, function(x){gIntersects(grid_pts_all, x, byid = TRUE)})
  keep <- apply(poss_C_pts, 1, function(x){all(x)})
  grid_pts <- SpatialPoints(grid_scale[keep,])
  n_poss_C_pts <- nrow(grid_pts@coords)
  
  #Compute areas in error (area_out)
  C_poss <- grid_pts@coords
  area_out <- matrix(nrow = n_poss_C_pts, ncol = n_train)
  for (i in 1:n_poss_C_pts) {
     l <- make_l(C = C_poss[i,], theta = thetas)
    for (j in 1:n_train) {
      pts_on_l_j <- pts_on_l(l = l, cont = train[[j]], under = FALSE)
      poss_j <- make_poly(pts_on_l_j, "test_cont")
      diff_reg1 <- gDifference(poss_j, train[[j]])
      diff_reg2 <- gDifference(train[[j]], poss_j)
      area_out[i, j] <- gArea(diff_reg1) + gArea(diff_reg2)
    }
  }
  
  #Find estimated center point
  opt_ind <- which.min(apply(area_out, 1, max))
  C_hat <- matrix(grid_pts@coords[opt_ind,], ncol = 2)
  
  #Make sets of lines, l, for each C_hat and get y's
  l_untrim <- make_l(C = C_hat, thetas)
  l <- lapply(l_untrim, function(x){gIntersection(x, box)})
  l_lengths <- sapply(l, function(x){as.numeric(gLength(x))})
  
  
  #priors
  mu0 <- (2/3)*l_lengths
  Lambda0_sigma2 <- .5^2
  Lambda0 <- Lambda0_sigma2*diag(p) 
  betaKappa0 <- 10
  betaSigma0 <- (l_lengths/2)/qnorm(.995)
  alphaAlpha0 <- 0
  betaAlpha0 <- 1
  
  #compute y's
  pts <- lapply(train, function(x) {pts_on_l(l, x, under = FALSE)})
  y <- sapply(pts, function(x){apply(x, 1, function(y){get_dist(y, C_hat)})})
  
  #initial values for MCMC
  mu_ini <- apply(y, 1, mean)
  kappa_ini <- 2#5
  sigma_obs <-  apply(y, 1, sd)
  sigma_ini <- sigma_obs
  sigma_ini[sigma_ini == 0] <- .0001
  alpha_ini <- 0.5
  
  #MCMC settings
  n_iter <- 30000 #40000
  burn_in <- 10000 #10000
  sigmaProp_sigma2 <- .05
  muProp_sigma2 <- .01
  kappaPropSD <- .01
  alphaPropSD <- .01
  cz <- contig_zero(sigma_obs)
  g_space <- 5
  g_start <- unlist(apply(cz, 1, function(x){seq(x[1], x[2], g_space)}))
  g_end <- c(g_start[2:length(g_start)] - 1, p)
  
  #non-scaler proposals for MCMC
  theta_dist_est <- theta_dist_mat(thetas)
  sigmaPropCov = sigmaProp_sigma2*compSigma(sigma_ini, kappa_ini, alpha_ini, theta_dist_est)
  muPropCov <- muProp_sigma2*compSigma(sigma_ini, kappa_ini, alpha_ini, theta_dist_est)
  
  #run MCMC
  start_time <- proc.time()
  fits <- ContouR::RunMCMC(nIter = n_iter, y,
                  mu = mu_ini, mu0 = mu0, Lambda0, muPropCov,
                  kappa = kappa_ini, betaKappa0, kappaPropSD,
                  sigma = sigma_ini, betaSigma0 = betaSigma0, sigmaPropCov,
                  gStart = g_start - 1, gEnd = g_end - 1,
                  thetaDist = theta_dist_est, 
                  alpha = alpha_ini, alphaPropSD, alphaAlpha0, betaAlpha0)
  end_time <- proc.time()
  
  #parameter estimates
  mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
  kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
  sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)
  alpha_est <- mean(fits$alpha[(burn_in + 1):n_iter])
  
  #check MCMC diagnostic
  # N_sigma <- N_mu <- rep(NA, p)
  # for (i in 1:p) {
  #   N_sigma[i] <- max(raftery.diag(fits$sigma[i,], q = .025, r = .0125)$resmatrix[,"N"],
  #                     raftery.diag(fits$sigma[i,], q = .975, r = .0125)$resmatrix[,"N"])
  #   N_mu[i] <- max(raftery.diag(fits$mu[i,], q = .025, r = .0125)$resmatrix[,"N"],
  #                   raftery.diag(fits$mu[i,], q = .975, r = .0125)$resmatrix[,"N"])
  # }
  # 
  # N_kappa <- max(raftery.diag(fits$kappa, q = .025, r = .0125)$resmatrix[,"N"],
  #                raftery.diag(fits$kappa, q = .975, r = .0125)$resmatrix[,"N"]) 
  
  
  #########################
  #Evaluate fits
  #########################
  # Sigma_obs <- cov(t(y))
  # Sigma_ini <- compSigma(sigma_ini, kappa_ini, alpha_ini, theta_dist_est)
  # Sigma_est <- compSigma(sigma_est, kappa_est,  alpha_ini, theta_dist_est)
  # par(mfrow = c(1, 2))
  # image.plot(Sigma_est, zlim = c(-.004, .013))
  # image.plot(Sigma_obs, zlim = c(-.004, .013))
  # 
  # 
  # plot(sigma_est, type= "l")
  # points(sigma_ini, type= "l", col = "blue")
  # points(betaSigma0, type = "l", lty = 2)
  # plot(mu_est, type= "l")
  # points(mu_ini, col= 'blue', type= "l") 
  # 
  # plot(fits$kappa, type= "l")
  # plot(fits$alpha, type= "l")
  
  
  #posterior field
  n_gen <- 100
  gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est, 
                    alpha = alpha_est, sigma = sigma_est, Cx = C_hat[1],
                    Cy = C_hat[2], thetas = thetas)
  prob <- prob_field(gens$polys, nrows = n_grid, ncols = n_grid)
  
  #compare results to generated contours
  pdf(sprintf("/users/hdirector/desktop/probAndGen%i.pdf", test_year),
      width = 8.5, height = 4)
  par(mfrow = c(1, 2))
  image.plot(prob, xaxt = "n", yaxt = "n", col = viridis(20),
             main = "Probability Field")
  for (i in train_ind) {
    points(coords_scale[[i]], col = 'black', type= "l")
  }
  points(coords_scale[[n_train + 1]],  col = 'red', type= "l")
  
  
  #creds <- cred_regs(prob, c(80, 90, 95), nrows = n_grid, ncols = n_grid)
  plot(creds[[3]], col = 'navy', lwd = 2, border = "navy", 
       main = "Credible Intervals")
  plot(creds[[2]], col = 'blue', lwd = 2, add = TRUE, border = "blue")
  plot(creds[[1]], col = 'lightblue', lwd = 2, add = TRUE, border = "lightblue")
  
  for (i in train_ind) {
    points(coords_scale[[i]],  type= "l")
  }
  points(coords_scale[[eval_ind]], type= "l", col= 'red')
  legend("topright",  fill = c("navy", "blue", "lightblue", "black", "red"),
         legend = c("80%","90%", "95", "Training Contours", "Test Contour"), 
         cex = .45)
  dev.off()
  
  #####################
  #formal evaluation
  #####################
  
  #Find C_hat for evaluation contour
  keep_eval <- gIntersects(grid_pts_all, eval, byid = TRUE)
  poss_pts_eval <- SpatialPoints(grid_scale[keep_eval,])
  n_poss_pts_eval <- nrow(grid_pts_eval@coords)
  C_poss_eval <- grid_pts_eval@coords
  
  area_out_eval <- rep(NA, n_test_pts_eval)
  for (i in 1:n_poss_pts_eval) {
    l <- make_l(C = C_poss_eval[i,], theta = thetas)
    pts_on_l <- pts_on_l(l = l, cont = eval, under = FALSE)
    poss <- make_poly(pts_on_l, "test_cont")
    diff_reg1 <- gDifference(poss, eval)
    diff_reg2 <- gDifference(eval, poss)
    area_out_eval[i] <- gArea(diff_reg1) + gArea(diff_reg2)
  }
  opt_ind_eval <- which.min(area_out_eval)
  C_hat_eval <- matrix(grid_pts_eval@coords[opt_ind_eval,], ncol = 2)
  
  
  #evaluate coverage
  pdf(sprintf("/users/hdirector/desktop/cred_ints%i.pdf", test_year),
              width = 8.5, height = 4)
  par(mfrow = c(1, 3))
  cover <- sapply(creds, 
         function(x){eval_cred_reg(truth = eval, cred_reg = x, center = C_hat_eval, 
                                   p_test = p_test, nrows = n_grid, ncols = n_grid,
                                   plotting = TRUE)})
  mtext(text = "80% Credible Interval", side = 3, line = -3, at = .18, outer = TRUE)
  mtext(text = "90% Credible Interval", side = 3, line = -3, at = .5, 
        outer = TRUE)
  mtext(text = "95% Credible Interval", side = 3, line = -3, at = .82,
        outer = TRUE)
  dev.off()
  
  save(cover, file = sprintf("/users/hdirector/desktop/cover_year%i.rda", test_year))
}