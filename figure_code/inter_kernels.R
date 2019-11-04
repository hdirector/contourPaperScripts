rm(list = ls())
library("ContouR")
n_obs <- 20
#TO DO (FIX PACKAGE, WHY IS IT LOADING THE WRONG DATA OBJECTS
load("data/stop_sign.rda")
load("data/tie.rda")
load("data/tree.rda")

#stop_sign
p <- length(tree$mu)
theta_space <- 2*pi/p
thetas <- seq(theta_space/2, 2*pi, by = theta_space)
obs_stop <- gen_conts(n_sim = n_obs, mu = stop_sign$mu, kappa = stop_sign$kappa,
                     sigma = stop_sign$sigma, Cx = stop_sign$Cx, Cy = stop_sign$Cy,
                     theta1 = thetas[1])
kern_stop <- find_inter_kernel(obs_stop$coords)

#tie
p <- length(tree$mu)
theta_space <- 2*pi/p
thetas <- seq(theta_space/2, 2*pi, by = theta_space)
n_obs <- 20
obs_tie <- gen_conts(n_sim = n_obs, mu = tie$mu, kappa = tie$kappa,
                      sigma = tie$sigma, Cx = tie$Cx, Cy = tie$Cy,
                      theta1 = thetas[1])
kern_tie <- find_inter_kernel(obs_tie$coords)


#tree
p <- length(tree$mu)
theta_space <- 2*pi/p
thetas <- seq(theta_space/2, 2*pi, by = theta_space)
obs_tree <- gen_conts(n_sim = n_obs, mu = tree$mu, kappa = tree$kappa,
                  sigma = tree$sigma, Cx = tree$Cx, Cy = tree$Cy,
                  theta1 = thetas[1])
kern_tree <- find_inter_kernel(obs_tree$coords)


pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/sample_shapes_kernels.pdf", 
    height = 4, width = 8.5)
par(mfrow = c(1, 3))
#stop_sign
plot(obs_stop$poly[[1]],xlim = c(.3, .7), ylim = c(.3, .7), lwd = .35,
     main = sprintf("%i Observed Contours and Shared Kernel", n_obs))
plot(kern_stop$poly, add = T, col = 'red', border = 'red')
for (i in 1:n_obs) {
  plot(obs_stop$polys[[i]], add = T, lwd = .35)
}

#tie
plot(obs_tie$poly[[1]],xlim = c(.3, .7), ylim = c(.3, .7), lwd = .35,
     main = sprintf("%i Observed Contours and Shared Kernel", n_obs))
plot(kern_tie$poly, add = T, col = 'red', border = 'red')
for (i in 1:n_obs) {
  plot(obs_tie$poly[[i]], add = T, lwd = .35)
}

#tree
plot(obs_tree$polys[[1]],xlim = c(.3, .7), ylim = c(.3, .7), lwd = .35,
     main = sprintf("%i Observed Contours and Shared Kernel", n_obs))
plot(kern_tree$poly, add = T, col = 'red', border  = 'red')
for (i in 1:n_obs) {
  plot(obs_tree$polys[[i]], add = T, lwd = .35)
}
dev.off()
