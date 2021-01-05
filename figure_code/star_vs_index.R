rm(list = ls())
set.seed(100)
library("ContouR")

#######################
#generate true contour
#######################
C <- matrix(data = c(.5, .5), ncol = 2)
p <- 50
theta <- sort(runif(50, 0, 2*pi))
cont_1 <- gen_conts(n_sim = 1, mu = shape2_2$mu, kappa = shape2_2$kappa,
                    sigma = shape2_2$sigma, C = shape2_2$C, 
                    thetas = theta)

########################################
#indexed sequence vs star-shaped indexed 
########################################
pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/star_v_point.pdf", 
    height = 2.5, width = 4)
par(mfrow = c(1, 2), oma = c(4, 0, 0, 0), mar = rep(0, 4))
plot(cont_1$polys[[1]], border = 'red', lwd = 2, col = 'lightblue')
points(cont_1$coords[[1]], pch = 20, cex = .5)

plot(cont_1$polys[[1]], border = 'red', lwd = 2, col = 'lightblue')
apply(cont_1$coords[[1]], 1, function(x){points(rbind(x, shape2_2$C), type= "l")})
points(C, pch = 3, col = 'purple', lwd = 3, cex = 2)
points(cont_1$coords[[1]], pch = 15, cex = .5, col = 'blue')

legend(-.5, .1, cex = .65, ncol = 2, xpd = NA, 
       legend = c(expression(paste("contour line (", bold(bar("S")), ")")),
                  expression(paste("point sequence (", bold("S"), ")")),
                  expression(paste("enclosed polygon (", bold(underline("S")), ")")),
                  expression(paste("starting point (", bold("C"), ")")),
                  expression(paste("lengths to contour (", bold(y), ")")),
                  expression(paste("star-shaped points (", bold(V), ")"))),
       pch = c(NA, 20, 15, 3, NA, 15), lty = c(1, NA, NA, 1, 1, 1, NA),
       col = c("red", "black", "lightblue", "purple", "black", "blue"),
       lwd = c(3, 3, NA, NA, 3, NA), pt.cex = 1, text.font = 1, fill = "white",
       bg = "white", border = "white")
dev.off()

###############################
#star-shaped approximation
################################
theta_space_2 <- 2*pi/p
theta_2 <- seq(theta_space_2/2, 2*pi, by = theta_space_2)
l <- make_l(C = C, theta = theta_2)
approx_pts <- pts_on_l(l = l, cont = cont_1$polys[[1]], under = TRUE)
approx <- make_poly(approx_pts, "approx_pts")

pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/approxContour.pdf", 
    height = 2.5, width = 4)
par(mfrow = c(1, 2), oma = c(4, 0, 0, 0), mar = rep(0, 4))
plot(cont_1$polys[[1]], border = 'red', lwd = 2, col = 'lightblue')
points(cont_1$coords[[1]], pch = 20, cex = .5)

plot(approx, border = 'red', lwd = 2, col = 'lightblue')
apply(approx_pts, 1, function(x){points(rbind(x, shape2_2$C), type= "l")})
points(C, pch = 3, col = 'purple', lwd = 3, cex = 2)
points(approx_pts, pch = 15, cex = .5, col = 'blue')

legend(-.5, .1, cex = .65, ncol = 2, xpd = NA, 
       legend = c(expression(paste("contour line (", bold(bar("S")), ")")),
                  expression(paste("point sequence (", bold("S"), ")")),
                  expression(paste("enclosed polygon (", bold(underline("S")), ")")),
                  expression(paste("starting point (", bold("C"), ")")),
                  expression(paste("lengths to contour (", bold(y), ")")),
                  expression(paste("star-shaped points (", bold(V), ")"))),
       pch = c(NA, 20, 15, 3, NA, 15), lty = c(1, NA, NA, 1, 1, 1, NA),
       col = c("red", "black", "lightblue", "purple", "black", "blue"),
       lwd = c(3, 3, NA, NA, 3, NA), pt.cex = 1, text.font = 1, fill = "white",
       bg = "white", border = "white")
dev.off()

