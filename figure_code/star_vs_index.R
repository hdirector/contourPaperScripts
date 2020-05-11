rm(list = ls())
set.seed(3)
library("ContouR")

#######################
#generate true contour
#######################
C <- matrix(data = c(.5, .5), ncol = 2)
p <- 75
theta_space <- 2*pi/p
theta <- seq(theta_space/2, 2*pi, by = theta_space)
cont <- gen_conts(n_sim = 1, mu = shape2_2$mu, kappa = shape2_2$kappa,
                  sigma = shape2_2$sigma, C = shape2_2$C, 
                  thetas = shape2_2$theta)
########################################
#indexed sequence vs star-shapedindexed 
########################################
pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/star_v_point.pdf", 
    height = 2.5, width = 4)
par(mfrow = c(1, 2), oma = c(4, 0, 0, 0), mar = rep(0, 4))
plot(cont$polys[[1]], border = 'red', lwd = 2, col = 'lightblue')
points(cont$coords[[1]], pch = 20, cex = .5)

plot(cont$polys[[1]], border = 'red', lwd = 2, col = 'lightblue')
apply(cont$coords[[1]], 1, function(x){points(rbind(x, shape2_2$C), type= "l")})
points(C, pch = 3, col = 'purple', lwd = 3, cex = 2)
points(cont$coords[[1]], pch = 15, cex = .5, col = 'blue')

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


#######################################
#Define generating lines, x, w, y, Cx, Cy
#######################################
#plot toy generating linesa
p <- 6; r <- .5
theta_space <- (2*pi)/p
thetas <- theta_space/2 + theta_space*(0:(p-1))
Cx <- .5; Cy <- .5
x_gen_line <- Cx + r*cos(thetas) 
y_gen_line <- Cy + r*sin(thetas)
#pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/method_components_demo.pdf")
plot(0, 0, xlim = c(-.05, 1.05), ylim = c(-.05, 1.05), col = 'white', xlab = "", ylab = "",
     xaxt = "n", yaxt = "n")
arrows(Cx, Cy, x_gen_line[1], y_gen_line[1])
for (i in 2:p) {
  arrows(Cx, Cy, x_gen_line[i], y_gen_line[i])
}
text(.97, .76, expression("\u2113"[1]))
text(.49, 1.03, expression("\u2113"[2]))
text(.043, .77, expression("\u2113"[3]))
text(.045, .24, expression("\u2113"[4]))
text(.5, -.03, expression("\u2113"[5]))
text(.97, .23, expression("\u2113"[6]))

#x, y, w
r_poly <- c(.2, .2, .35, .3, .2, .1)
x <- cbind(Cx + r_poly*cos(thetas), Cy + r_poly*sin(thetas))
#plot poly
points(rbind(x, x[1,]), type = "l", col = 'blue', lwd = 2)


#x and z
text(.72, .59, labels = expression(x[1]))
points(x[1,1], x[1,2], col = 'purple', pch = 20)
points(rbind(c(Cx, Cy), x[1,]), type = "l", col = 'purple', 
       lwd = 2)
text(.57, .58, expression(z[1]))


#plot Cx, Cy
points(Cx, Cy, col = 'red', pch = 20, cex = 1.5)
text(Cx - .02, Cy - .06, "C")
#dev.off()

# 
# #####################
# #Show how to find C
# #####################
# #pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/find_C.pdf",
# #    height = 4, width = 10)
# r1_poly <- c(.26, .2, .35, .26, .2, .12)
# r2_poly <- c(.25, .3, .3, .35, .15, .27)
# r3_poly <- c(.1, .1, .25, .3, .24, .22)
# r4_poly <- c(.2, .2, .2, .32, .1, .05)
# x1 <- cbind(Cx + r1_poly*cos(thetas), Cy + r1_poly*sin(thetas))
# x2 <- cbind(Cx + r2_poly*cos(thetas), Cy + r2_poly*sin(thetas))
# x3 <- cbind(Cx + r3_poly*cos(thetas), Cy + r3_poly*sin(thetas))
# x4 <- cbind(Cx + r4_poly*cos(thetas), Cy + r4_poly*sin(thetas))
# 
# #two contours
# par(mfrow = c(1, 3))
# plot(0, 0, xlim = c(.15, .75), ylim = c(.25, .8), col = 'white', xlab = "", ylab = "",
#      xaxt = "n", yaxt = "n")
# #lines
# points(rbind(x1, x1[1,]), type = "l", col = 'blue', lwd = 2)
# points(rbind(x2, x2[1,]), type= "l", col = 'pink', lwd = 2)
# 
# plot(0, 0, xlim = c(.15, .75), ylim = c(.25, .8), col = 'white', xlab = "", ylab = "",
#      xaxt = "n", yaxt = "n")
# #lines
# points(rbind(x1, x1[1,]), type = "l", col = 'blue', lwd = 2)
# points(rbind(x2, x2[1,]), type= "l", col = 'pink', lwd = 2)
# 
# #points
# points(rbind(x1, x1[1,]), pch = 20, col = 'blue')
# points(rbind(x2, x2[1,]), pch = 20, col = 'pink')
# 
# 
# line_eq <- function(p1, p2) {
#   m <- (p2[2] - p1[2])/(p2[1] - p1[1])
#   b <- p1[2] - m*p1[1]
#   return(list("m" = m, "b" = b))
# }
# 
# #connected line 1
# l1 <- line_eq(x1[1,], x2[1,])
# x_line1 <- seq(.2, x1[1,1], length = 10)
# y_line1 <- l1$m*x_line1 + l1$b
# points(x_line1, y_line1, type = "l", lwd = 2, lty = 2)
# 
# #connected line 2 (happens to be vertical)
# points(cbind(rep(.5, 2), c(.8, .26)), type = "l", lwd = 2, lty = 2)
# 
# #connected line 3
# l3 <- line_eq(x1[3,], x2[3,])
# x_line3 <- seq(.196, .73, length = 10)
# y_line3 <- l3$m*x_line3 + l3$b
# points(cbind(x_line3, y_line3), type = "l", lwd = 2, lty = 2)
# 
# points(Cx, Cy, col = 'red', pch = 20, cex = 1.5)
# text(Cx - .015, Cy - .025, "C")
# 
# 
# 
# #add more contours
# plot(0, 0, xlim = c(.15, .75), ylim = c(.25, .8), col = 'white', xlab = "", ylab = "",
#      xaxt = "n", yaxt = "n")
# 
# #lines
# points(rbind(x1, x1[1,]), type = "l", col = 'blue', lwd = 2)
# points(rbind(x2, x2[1,]), type= "l", col = 'pink', lwd = 2)
# points(rbind(x3, x3[1,]), type= "l", col = 'green', lwd = 2)
# points(rbind(x4, x4[1,]), type= "l", col = 'purple', lwd = 2)
# #points
# points(rbind(x1, x1[1,]), pch = 20, col = 'blue')
# points(rbind(x2, x2[1,]), pch = 20, col = 'pink')
# points(rbind(x3, x3[1,]), pch = 20, col = 'green')
# points(rbind(x4, x4[1,]), pch = 20, col = 'purple')
# 
# #Add conntected lines
# points(x_line1, y_line1, type = "l", lwd = 2, lty = 2)
# points(cbind(rep(.5, 2), c(.8, .26)), type = "l", lwd = 2, lty = 2)
# points(cbind(x_line3, y_line3), type = "l", lwd = 2, lty = 2)
# 
# points(Cx, Cy, col = 'red', pch = 20, cex = 1.5)
# text(Cx - .015, Cy - .025, "C")
# dev.off()
# 
# 
# ###########################
# #Can you have more than C?
# ##########################
# #pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/multiple_Cs.pdf",
# #    height = 4, width = 7)
# par(mfrow = c(1, 3))
# plot(0, 0, xlim = c(.15, .75), ylim = c(.25, .8), col = 'white', xlab = "", ylab = "",
#      xaxt = "n", yaxt = "n")
# points(rbind(x1, x1[1,]), type = "l", col = 'blue', lwd = 2)
# points(rbind(x1, x1[1,]), col = 'blue', pch = 20, cex  = 2)
# arrows(Cx, Cy, x1[1, 1], x1[1, 2])
# for (i in 2:p) {
#   arrows(Cx, Cy, x1[i, 1], x1[i, 2])
# }
# points(Cx, Cy, col = 'red', pch = 20, cex = 2)
# 
# Cx2 <- .4; Cy2 <- .6
# plot(0, 0, xlim = c(.15, .75), ylim = c(.25, .8), col = 'white', xlab = "", ylab = "",
#      xaxt = "n", yaxt = "n")
# points(rbind(x1, x1[1,]), type = "l", col = 'blue', lwd = 2)
# points(rbind(x1, x1[1,]), col = 'blue', pch = 20, cex  = 2)
# arrows(Cx2, Cy2, x1[1, 1], x1[1, 2])
# for (i in 2:p) {
#   arrows(Cx2, Cy2, x1[i, 1], x1[i, 2])
# }
# points(Cx2, Cy2, col = 'red', pch = 20, cex = 2)
# 
# 
# Cx3 <- .4; Cy3 <- .4
# plot(0, 0, xlim = c(.15, .75), ylim = c(.25, .8), col = 'white', xlab = "", ylab = "",
#      xaxt = "n", yaxt = "n")
# points(rbind(x1, x1[1,]), type = "l", col = 'blue', lwd = 2)
# points(rbind(x1, x1[1,]), col = 'blue', pch = 20, cex  = 2)
# arrows(Cx3, Cy3, x1[1, 1], x1[1, 2])
# for (i in 2:p) {
#   arrows(Cx3, Cy3, x1[i, 1], x1[i, 2])
# }
# points(Cx3, Cy3, col = 'red', pch = 20, cex = 2)
# #dev.off()
# 
# 
# 
# 
