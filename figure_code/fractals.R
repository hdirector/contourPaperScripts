library("pracma")
a <- fractalcurve(8, "snowflake")
n_plot <- length(a$x)

circle_pts <- function(x, y, delta = .022, N = 100) {
  theta <- seq(0, 2*pi, length = N)
  return(cbind(x + delta*cos(theta), y + delta*sin(theta)))
}

pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Figures/fractals.pdf",
    height = 3, width = 3)
layout(matrix(nrow = 5, ncol = 1, data = c(1, 1, 1, 2, 2)))
par(oma = rep(0, 4), mar = rep(0, 4))

#fractal with circle
plot(cbind(a$x, a$y), type= 'l', xaxt = "n", yaxt = "n",
     xlab = "", ylab = '', xlim = c(-.05, 1.05), ylim = c(-.95, .3),
     bty = 'n')

space <- 5
centers <- matrix(nrow = 0, ncol = 2)
for (i in seq(space/2, n_plot, space)) {
    x_temp <- mean(a$x[i:(i + 1)])
    y_temp <- mean(a$y[i:(i + 1)])
    points(circle_pts(x_temp, y_temp), type= 'l', col = 'darkgrey')
    centers <- rbind(centers, c(x_temp, y_temp))
}
points(centers, col = 'red', type = 'l', lwd = 2)
points(cbind(a$x, a$y), type= 'l',lwd = 1.2)

#zoomed in
n_zoom <- 64
plot(cbind(a$x[1:n_zoom], a$y[1:n_zoom]), type= 'l', xaxt = "n", yaxt = "n",
     xlab = "", ylab = '', xlim = c(0, 1/3), ylim = c(-.14, .2),
     bty = 'n')
space <- 4
centers_zoom <- matrix(nrow = 0, ncol = 2)
for (i in seq(space/2, n_zoom, space)) {
  x_temp <- mean(a$x[i:(i + 1)])
  y_temp <- mean(a$y[i:(i + 1)])
  points(circle_pts(x_temp, y_temp), type= 'l', col = 'darkgrey', 
         lwd = 2)
  centers_zoom <- rbind(centers_zoom, c(x_temp, y_temp))
}
points(cbind(a$x[1:n_zoom], a$y[1:n_zoom]), type= 'l', lwd = 2)
points(centers_zoom, col = 'red', type = 'l', lwd = 2)
points(centers_zoom, col = 'blue', pch = 15, cex = .8)
legend("bottom", legend = c("fractal contour", "covering circles",
                            "circle centers", "piecewise linear line"),
       col = c("black", "darkgrey", "blue", "red"),
       lty = c(1, 1, NA, 1), ncol = 2, cex = .75,
       text.font = 2, lwd = 2, bty = 'n', pch = c(NA, NA, 15, NA))
dev.off()

