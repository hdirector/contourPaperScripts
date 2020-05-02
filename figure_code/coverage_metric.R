rm(list = ls())
library("ContouR")
set.seed(3)
n_grid <- 150

#######################
#generate true contour
#######################
C <- matrix(data = c(.5, .5), ncol = 2)
p <- 50
theta_space <- 2*pi/p
theta <- seq(theta_space/2, 2*pi, by = theta_space)
cont <- gen_conts(n_sim = 1, mu = shape2_2$mu, kappa = shape2_2$kappa,
                  sigma = shape2_2$sigma, C = shape2_2$C, 
                  thetas = shape2_2$theta)

#############################
#make up a credible interval
#############################
outer <- gen_conts(n_sim = 1, mu = shape2_2$mu + .12 + rnorm(p, 0, .005) , 
                   kappa = shape2_2$kappa, sigma = shape2_2$sigma, 
                   C = shape2_2$C,  thetas = shape2_2$theta)
inner <- gen_conts(n_sim = 1, mu = shape2_2$mu - .07 + rnorm(p, 0, .005), 
                   kappa = shape2_2$kappa, sigma = shape2_2$sigma,
                   C = shape2_2$C,  thetas = shape2_2$theta)
cred <- gDifference(outer$polys[[1]], inner$polys[[1]])

###########################
#make figure
###########################
pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/metric_demo.pdf", 
    height = 6.1, width = 8)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
p_test <- 100
theta_space <- 2*pi/p_test
thetas_evals <- seq(theta_space, 2*pi, theta_space)
eval_cred_reg(truth = cont$polys[[1]], cred_reg = cred, center = C,
              thetas = thetas_evals,nrows = n_grid, ncols = n_grid,
              plotting = TRUE)
legend("bottom",   pch = c(NA, 15, 3, NA, NA), lty = c(1, NA, NA, 1, 1),
       lwd = c(5, NA, 5, 5, 5), col = c("red", "lightcyan2", "darkgreen", "black", 
                                        "blue"),
       ncol = 2, pt.cex = c(1, 2, 1,1,1),
       legend = c(expression(paste("contour line (", bold(bar("S"))[i], ")")),
                  expression(paste("credible interval (", bold(I [1 - alpha]), ")")), 
                  expression(paste("starting point (", bold("C"),  ")")),
                  expression(paste("intersection of ", I [1 - alpha], 
                                   " and line k (",
                                   W["i,k"], " = 1)")),
                  expression(paste("intersection of ", I [1 - alpha], 
                                   " and line k (",
                                   W["i,k"], " = 0)"))))
dev.off()

