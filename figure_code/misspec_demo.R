#Assess coverage under particular conditions
rm(list = ls())
set.seed(103)
library("ContouR")
library("coda")

#Read in simulation settings
gens <- list()
for (task_id in 1:3) {
  task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/misspecTaskTable.rda"
  load(task_path)
  task <- task_table[task_id,]
  print(sprintf("task_id: %i", task_id))
  attach(task)
  
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
  
  #generate two contours for demonstrations
  n_obs <- 1
  box <- bbox()
  bd <- box@polygons[[1]]@Polygons[[1]]@coords
  gens[[task_id]] <- gen_misspec(n_sim = n_obs, mu = mu_true, kappa = kappa_true,
                      sigma = sigma_true, C = C_true, thetas = thetas_true, 
                      r1_min = r1_min, r1_max = r1_max, r2_min = r2_min,
                      r2_max = r2_max, n_curl_min = n_curl_min,
                      n_curl_max = n_curl_max, bd = bd)
}
pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/misspec_exs.pdf",
    height = 2, width = 4.5)
par(mfrow = c(1, 3), oma = rep(0, 4), mar = rep(0, 4))
plot(gens[[1]]$polys[[1]], xlim = c(0.1, 0.9), ylim = c(0.1, 0.9))
plot(gens[[2]]$polys[[1]], xlim = c(0.1, 0.9), ylim = c(0.1, 0.9))
plot(gens[[3]]$polys[[1]], xlim = c(0.1, 0.9), ylim = c(0.1, 0.9))
dev.off()
