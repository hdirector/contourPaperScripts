#generate a star-shaped polygon

rm(list = ls())
library("sp")
library("MASS")
library("viridis")
library("rgeos")
source('R/planes_intersections_polys.R')
source('R/misc.R')

#true circular contour
p <- 20
theta <- seq(0, 2*pi, length = p + 1)
theta <- theta[1:p] #don't want repeated points 0/pi
mu <- c(seq(.1, .3, length = p/2), seq(.3, .1, length = p/2))
#mu <- rep(.05, p)
center <- c(.5, .5)
truth <- cbind(center[1] + mu*cos(theta), center[2] + mu*sin(theta))
Sigma_sqExp <- .001*exp(-dist_mat_circle(p)) #true covariance


#Generate contours
n_obs <- 200
y_obs <-  mvrnorm(n_obs, mu, Sigma_sqExp)
cut_off <- .05
y_obs[y_obs < cut_off] <- cut_off
  
#view observations and truth
plot(rbind(truth, truth[1,]), xlim = c(0, 1), ylim = c(0,1), type= "l", col = 'green')
obs_coords <- list()
for (i in 1:n_obs) {
  obs_coords[[i]] <- cbind(center[1] + y_obs[i,]*cos(theta), 
                           center[2] + y_obs[i,]*sin(theta))
  points(rbind(obs_coords[[i]], obs_coords[[i]][1,]), type= "l")
}


bb <- bbox()
print_nums <- c(1, 10, 25, 50, 100, 200)
#pdf("Figures/intersection_kernel.pdf")
par(mfrow = c(2, 3))
for (i in 1:n_obs) {
  kernel_i <- find_kernel(obs_coords[[i]])
  if (i == 1) {
    shared_kernel <- kernel_i
  } else {
    shared_kernel <- gIntersection(shared_kernel, kernel_i)
  }
  stopifnot(!is.null(kernel_i))
  if (i %in% print_nums) {
    plot(bb,  main = sprintf("%i Contours", i))
    plot(shared_kernel, add = T, col = 'blue', border = "blue")
    for (j in 1:i) {
      points(rbind(obs_coords[[j]], obs_coords[[j]][1,]), type= "l")
    }
  }
}
#dev.off()
plot(bb)
plot(shared_kernel, add = T, col = 'red')
points(center[1], center[2], col = 'green')
plot(gCentroid(shared_kernel), add = T, col = 'green')



