rm(list = ls())
library("pracma")
library("ContouR")

#make fractal in [.1, .9] x [.1, .9] box
frac_orig <- fractalcurve(8, "snowflake")
frac_orig_list <- list()
frac_orig_list[[1]] <- cbind(frac_orig$x, frac_orig$y)
frac_rescale <- rescale(coords = frac_orig_list, eps = .1, box_size = .01)
frac <- frac_rescale$coords_scale[[1]]

n_pts <- nrow(frac)
circle_pts <- function(x, y, delta, N = 100) {
  theta <- seq(0, 2*pi, length = N)
  return(cbind(x + delta*cos(theta), y + delta*sin(theta)))
}

pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Figures/fractals.pdf",
    height = 4, width = 3)
layout(matrix(nrow = 7, ncol = 1, data = c(1, 1, 1, 2, 2, 2, 3)))
par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))

#fractal with circle
plot(frac, type= 'l', xaxt = "n", yaxt = "n",
     xlab = "", ylab = '', xlim = c(0, 1), ylim = c(0, 1),
     bty = 'n')

space <- 5
delta <- .015
centers <- matrix(nrow = 0, ncol = 2)
secs <- seq(ceiling(space/2), n_pts, space)
n_secs <- length(secs)
for (i in secs) {
  if (i != secs[n_secs]) {
    ind <- i:(i + space - 1)
  } else { #edge case at end
    loop_add <- space - length(secs[n_secs]:n_pts)
    ind <- c(secs[n_secs]:n_pts, 1:loop_add)
  }
  stopifnot(length(ind) == space)
  x_temp <- mean(frac[ind, 1])
  y_temp <- mean(frac[ind, 2])
  points(circle_pts(x_temp, y_temp, delta = delta), type= 'l', 
           col = 'darkgrey')
  centers <- rbind(centers, c(x_temp, y_temp))
}
n_centers <- nrow(centers)
points(rbind(centers, centers[1,]), col = 'red', type = 'l', lwd = 1)
points(frac, type= 'l', lwd = 1)

#compute star-shaped approximation to centers
frac_cont <- list()
frac_cont[[1]] <- make_poly(coords = frac, name = "fractal")
bd <- bbox(eps = 0)
p <- 200
theta_space <- 2*pi/p
thetas <- seq(theta_space/2, 2*pi, length = p + 1)
centers_cont <- list()
centers_cont[[1]] <- make_poly(centers, "centers")
C_est <- best_C(bd = bd, conts = frac_cont, thetas = thetas, space = .04,
                parallel = TRUE)
centers_star <- error_areas(conts = centers_cont, C = C_est, thetas = thetas,
                            under = FALSE, ret_all = TRUE)
plot(centers_star$approxs[[1]], add = T, border = 'purple', lwd = 1)

diff_area <- gArea(gDifference(centers_star$approxs[[1]], frac_cont[[1]])) + 
  gArea(gDifference(frac_cont[[1]], centers_star$approxs[[1]]))
diff_area/gArea(frac_cont[[1]])

#zoomed in
n_zoom <- 64
plot(frac[1:n_zoom,], type= 'l', xaxt = "n", yaxt = "n",
     xlab = "", ylab = '', xlim = c(.07, .37), ylim = c(.59, .89),
     bty = 'n')
cent_ind <- c(n_centers, 1:(n_zoom/space + 1))
points(centers[cent_ind,], col = 'red', type = 'l', lwd = 1)
points(centers[cent_ind,], col = 'blue', pch = 15, cex = .8)
for (i in cent_ind) {
  points(circle_pts(centers[i,1], centers[i,2], delta = delta), type= 'l', 
         col = 'darkgrey')
}
points(centers_star$approxs[[1]]@polygons[[1]]@Polygons[[1]]@coords[120:135,],
       col = 'purple', type= 'l')
points(frac[1:n_zoom,], type= 'l')
plot(0,0, col = 'white', bty = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")
legend("center", legend = c(expression(paste("fractal contour (", bold("F"), ")")), 
                            expression(paste(bold(delta), "-cover ({U})")),
                            expression(paste("circle centers ({", C[i], "})")),
                            expression(paste("contour line (", bold(bar("S")), ")")),
                            "star-shaped representation"),
       col = c("black", "darkgrey", "blue", "red", "purple"),
       lty = c(1, 1, NA, 1, 1), ncol = 2, cex = .85,
       text.font = 2, lwd = 2, bty = 'n', pch = c(NA, NA, 15, NA, NA))
dev.off()


