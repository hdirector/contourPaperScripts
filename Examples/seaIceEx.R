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
start_year <- 2008
end_year <- 2017
xmn = -3850; xmx = 3750; ymn = -5350; ymx = 5850
box_size <- 25

#setting for executing 
eps <- .1 
n_gen <- 100
p_init <- 110
p_eval <- 100
delta <- .05
space <- .025
step <- .1
under <- FALSE

#MCMC settings
n_iter <- 50000
burn_in <-  25000
sigmaProp_sigma2 <- .05
muProp_sigma2 <- .05
kappaProp_SD <- .3

#derived settings
years <- start_year:end_year
n_years <- length(years)
n_train <- n_years - 1
box <- bbox(eps) 

#make a grid of points at the center of each square 
x_cent <- seq(xmn + box_size/2, xmx, box_size)
y_cent <- seq(ymn + box_size/2, ymx, box_size)
grid <- expand.grid(x_cent, y_cent)

#regions in and around central Arctic
reg <- reg_info$regions[[1]]
reg <- rm_holes(raster::aggregate(reg))
reg_coords <- reg@polygons[[1]]@Polygons[[1]]@coords

#load September data
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
                        land = land, grid = grid, bd = reg_coords)

coords_scale <- temp_rescale$coords_scale
grid_scale <- temp_rescale$grid_scale
bd_scale <- temp_rescale$bd_scale
bd_scale <- bd_scale[nrow(bd_scale):1, ] #reverse order to match generated contours
box_size_scale <- temp_rescale$box_size

#make observed contour polgons
obs <- list()
obs$polys <- lapply(coords_scale, function(x){make_poly(x, "contours")})

#make shifted land polygons
box <- bbox(0)
land_scale <- temp_rescale$land_scale
land_scale <- lapply(land_scale, function(x){make_poly(x, "land")})
land_scale <- lapply(land_scale, function(x){keep_poly(gIntersection(box, x))})
land_scale <- land_scale[sapply(land_scale, function(x){!is.null(x)})]

#make background water object
not_reg <- gDifference(box, make_poly(bd_scale, "bd"))
for (i in 1:length(land_scale)) {
  not_reg <- gDifference(not_reg, land_scale[[i]])
}

#identify section with no variability  
bd_poly <- make_poly(bd_scale, "bd_poly")
no_var <- gIntersection(as(bd_poly, "SpatialLines"), 
                        as(obs$polys[[1]], "SpatialLines"))
for (i in 2:n_years) {
  no_var <- gIntersection(no_var, obs$polys[[i]])
}

#Find C and thetas for evaluation contour
theta_space_eval <- 2*pi/p_eval
thetas_eval_all <- seq(theta_space_eval/2, 2*pi, by = theta_space_eval)
C_eval <- best_C(bd = bd_scale, conts = obs$polys, thetas = thetas_eval_all,
                 space = space, parallel = TRUE, under = FALSE)
area_err <- mean(error_areas(conts = obs$polys, C = C_eval, thetas = thetas_eval_all, 
                under = FALSE))

l_untrim_eval <- make_l(C = C_eval, thetas_eval_all)
l_eval <- lapply(l_untrim_eval, function(x){gIntersection(x, bd_poly)})
on_can_arch_eval <- sapply(l_eval, function(x){gIntersects(no_var, x)})
l_eval <- l_eval[!on_can_arch_eval]
thetas_eval <- thetas_eval_all[!on_can_arch_eval]

#input settings
for (test_year in start_year:end_year)  {
  #Seperate out training and testing set
  eval_ind <- which(years == test_year)
  train <- obs
  train$polys[[eval_ind]] <- NULL
  eval <- list()
  eval$polys <- obs$polys[[eval_ind]]
  
  #find estimated center point, p, and theta
  ests <- find_CP(conts = train, delta, p_init, space, step, misspec = TRUE, 
                  parallel = TRUE)
  C_est <- ests$C_est
  thetas_est_all <- ests$thetas_est
  p_est_all <- length(ests$thetas_est)

  #Make sets of lines, l, removing no variability sections on Can. Arch.
  l_untrim <- make_l(C = C_est, thetas_est_all)
  l <- lapply(l_untrim, function(x){gIntersection(x, bd_poly)})
  on_can_arch <- sapply(l, function(x){gIntersects(no_var, x)})
  l_lengths_all <- sapply(l, function(x){as.numeric(gLength(x))})
  l_rm <- which(on_can_arch)
  n_rm <- length(l_rm)
  l_keep <- which(!on_can_arch)
  l <- l[l_keep]
  
  #adjust l, theta, p to account for Can. Arch.
  l_lengths <- l_lengths_all[l_keep]
  thetas_est <- thetas_est_all[l_keep]
  p_est <- p_est_all - n_rm
  
  #priors
  mu0 <- rep(.5*mean(l_lengths), p_est)
  Lambda0_sigma2 <- .5^2
  Lambda0 <- Lambda0_sigma2*diag(p_est) 
  betaKappa0 <- 10
  betaSigma0 <- (l_lengths/2)/qnorm(.995)
  
  #compute y's
  pts <- lapply(train$polys, function(x) {pts_on_l(l, x, under = FALSE)})
  y <- sapply(pts, function(x){apply(x, 1, function(y){get_dist(y, C_est)})})
  
  #initial values for MCMC
  mu_ini <- apply(y, 1, mean)
  kappa_ini <- .5
  sigma_obs <-  apply(y, 1, sd)
  sigma_ini <- sigma_obs
  sigma_ini[sigma_ini == 0] <- .0001

  #non-scaler proposals for MCMC
  theta_dist_est <- theta_dist_mat(thetas_est)
  sigmaPropSD = sqrt(diag(sigmaProp_sigma2*compSigma(sigma_ini, kappa_ini,
                                                     theta_dist_est)))
  muPropSD <- sqrt(diag(muProp_sigma2*compSigma(sigma_ini, kappa_ini, 
                                                theta_dist_est)))
  #run MCMC
  start_time <- proc.time()
  fits <- ContouR::RunMCMC(nIter = n_iter, y = y, mu = mu_ini, mu0 = mu0,
                           Lambda0 = Lambda0, muPropSD = muPropSD,
                           kappa = kappa_ini, betaKappa0, kappaProp_SD,
                           sigma = sigma_ini, betaSigma0 = betaSigma0, 
                           sigmaPropSD = sigmaPropSD, thetaDist = theta_dist_est)
  end_time <- proc.time()
  run_time <- end_time - start_time
  print(sprintf("fitted year %i complete, time: %s", test_year, run_time[3]))
  save(fits, file = sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Examples/ex_results/fits_year%i.rda", 
                            test_year))

  #parameter estimates
  mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
  kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
  sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)

  #posterior field
  gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est, 
                    sigma = sigma_est, C = C_est, thetas = thetas_est_all,
                    bd = bd_scale, fix_ind = l_rm, rand_ind = l_keep,
                    fix_y =  l_lengths_all[l_rm])

  #set up grid to align with data
  grid_need <- 1/box_size_scale
  xmn <- ymn <- 0 
  xmx <- ymx <- 1 
  x_bds <- seq(xmn, xmx, by = box_size_scale[1])
  y_bds <- seq(ymn, ymx, by = box_size_scale[2])

  #compute probabilities and credible intervals
  n_grid_y <- length(y_bds) - 1
  n_grid_x <- length(x_bds) - 1
  prob <- prob_field(polys = gens$polys, nrows = n_grid_y, ncols = n_grid_x)
  creds <- cred_regs(prob, cred_eval = c(80, 90, 95), nrows = n_grid_y,
                     ncols = n_grid_x)
  
  #####################
  #formal evaluation
  #####################
  #evaluate coverage
  cover <- sapply(creds, 
         function(x){eval_cred_reg(truth = eval$polys, cred_reg = x, 
                                   center = C_eval,  thetas = thetas_eval, 
                                   nrows = n_grid, ncols = n_grid)})
  
  if (test_year == 2017) {
    pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/seaIce_cred.pdf",
        height = 3.6, width = 3.7)
    par(oma = rep(0, 4), mar = c(2, 0, 0, 0))
    temp <- eval_cred_reg(truth = eval$polys, cred_reg = creds[[1]], 
                          center = C_eval,  thetas = thetas_eval, 
                          nrows = n_grid, ncols = n_grid, plotting = TRUE,
                          land = land_scale, not_reg  = not_reg)
    legend(-.07, -.01, cex = .48, ncol = 2,
           legend = c(expression(paste("contour line (", bold(bar("S"))[i], ")")), 
                      expression(paste("credible interval (", bold(I [1 - alpha]), ")")),
                      "land", "outside region",
                      expression(paste("starting point (", bold("C"),  ")")),
                      expression(paste("intersection of ", bold(I [1 - alpha]), 
                                       " and line k (",
                                       W["i,k"], " = 1)")),
                      expression(paste("intersection of ", bold(I [1 - alpha]), 
                                       " and line k (",
                                       W["i,k"], " = 0)"))),
           pch = c(NA, rep(22, 3), 3, NA, NA), lty = c(1, rep(NA, 4), 1, 1),
           col = c("red", rep(NA, 3), "darkgreen", "black", "blue"),
           lwd = c(3, rep(NA, 2), 3, 3, 3), pt.cex = 1,
           text.font = 1, fill = c(NA, "lightcyan2",  "grey", "beige", NA,  NA, NA),
           xpd = NA,  bg = "white", border = c("white", rep("black",3),"white", "white"))
    dev.off()
  }
 
  
  save(cover, file = sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Examples/ex_results/cover_year%i.rda", 
                             test_year))
  print(sprintf("completed analysis of year %i", test_year))
}
