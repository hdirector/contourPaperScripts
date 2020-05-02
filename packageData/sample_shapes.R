rm(list = ls())
library("viridis")
library("fields")
library("ContouR")
set.seed(103)

##############################
#constants across shapes
##############################
n_grid <- 200
n_gen <- 100
p <- 50
C <- c(.5, .5)
theta_space <- 2*pi/p
theta <- seq(theta_space/2, 2*pi, by = theta_space)

############################
#shape A
############################
#coordinates and polygon of mean
r <- c(seq(.29, .35, length = 4), seq(.35, .29, length = 4),
       seq(.29, .38, length = 4), seq(.38, .25, length = 4),
       seq(.25, .17, length = 5), seq(.17, .25, length = 5),
       rep(.25, 10), seq(.25, .15, length = 5), seq(.15, .25, length = 5), 
       rep(.25, 4))
template_shapeA <- cbind(C[1] + r*cos(theta), C[2] + r*sin(theta))
template_shapeA <- rescale(coords = list(template_shapeA), eps = .2, box_size = 1) #eps = .2

#shape parameters
shapeA_poly <- make_poly(template_shapeA$coords_scale, "shape_2")

l_untrim <- make_l(C = C, theta)
l <- lapply(l_untrim, function(x){gIntersection(x, shapeA_poly)})
muA <- sapply(l, function(x){as.numeric(gLength(x))})
kappaA_1 <- 1
kappaA_2 <- 2
kappaA_4 <- 4
sigmaA <- c(seq(.035, .08, length = floor(p/8)), 
            seq(.08, .035, length = ceiling(p/8)), rep(.035, p/4),
            seq(.035, .08, length = floor(p/8)), 
            seq(.08, .035, length = ceiling(p/8)), rep(.035, p/4))


#generate probability distributions
gens_shapeA_1 <- gen_conts(n_sim = n_gen, mu = muA, kappa = kappaA_1,
                         sigma = sigmaA, C, thetas = theta)
probA_1 <- prob_field(polys = gens_shapeA_1$polys, nrows = n_grid, ncols = n_grid)
gens_shapeA_2 <- gen_conts(n_sim = n_gen, mu = muA, kappa = kappaA_2,
                           sigma = sigmaA, C, thetas = theta)
probA_2 <- prob_field(polys = gens_shapeA_2$polys, nrows = n_grid, ncols = n_grid)
gens_shapeA_4 <- gen_conts(n_sim = n_gen, mu = muA, kappa = kappaA_4,
                           sigma = sigmaA, C, thetas = theta)
probA_4 <- prob_field(polys = gens_shapeA_4$polys, nrows = n_grid, ncols = n_grid)

#save shape parameters
shapeA_1 <- list("mu" = muA, "kappa" = kappaA_1, "sigma" = sigmaA, "C" = C, 
               "theta" = theta)
shapeA_2 <- list("mu" = muA, "kappa" = kappaA_2, "sigma" = sigmaA, "C" = C, 
                 "theta" = theta)
shapeA_4 <- list("mu" = muA, "kappa" = kappaA_4, "sigma" = sigmaA, "C" = C, 
                 "theta" = theta)
save(shapeA_1, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/shapeA_1.rda")
save(shapeA_2, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/shapeA_2.rda")
save(shapeA_4, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/shapeA_4.rda")


#########################
#shape B
########################
#coordinates and polygon of mean
r <- .3
template_shapeB <- cbind(C[1] + r*cos(theta), C[2] + r*sin(theta))
shapeB_poly <- make_poly(template_shapeB, "shape_1")

#shape parameters
l_untrim <- make_l(C = C, theta)
l <- lapply(l_untrim, function(x){gIntersection(x, shapeB_poly)})
muB <- sapply(l, function(x){as.numeric(gLength(x))})
kappaB <- 2
sigmaB <-  sigmaA


#generate probability distribution
gens_shapeB <- gen_conts(n_sim = n_gen, mu = muB, kappa = kappaB,
                         sigma = sigmaB, C, thetas = theta)
probB <- prob_field(polys = gens_shapeB$polys, nrows = n_grid, ncols = n_grid)

#save shape parameters
shapeB <- list("mu" = muB, "kappa" = kappaB, "sigma" = sigmaB, "C" = C, 
               "theta" = theta)
save(shapeB, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/shapeB.rda")


##########################
#shape C 
##########################
#coordinates and polygon of mean
r <- c(seq(.3, .15, length = 5), seq(.15, .3, length = 5),
       seq(.29, .15, length = 4), seq(.15, .3, length = 4),
       seq(.3, .17, length = 6), seq(.17, .32, length = 6),
       seq(.32, .23, length = 5), seq(.15, .35, length = 5),
       seq(.36, .17, length = 5), seq(.17, .3, length = 5))
template_shapeC <- cbind(C[1] + r*cos(theta), C[2] + r*sin(theta))

#shape parameters
shapeC_poly <- make_poly(template_shapeC, "shape_3")
l_untrim <- make_l(C = C, theta)
l <- lapply(l_untrim, function(x){gIntersection(x, shapeC_poly)})
muC <- sapply(l, function(x){as.numeric(gLength(x))})
kappa_C <- 2
sigmaC <- sigmaA

#generate probability distribution
gens_shapeC <- gen_conts(n_sim = n_gen, mu = muC, kappa = kappa_C,
                         sigma = sigmaC, C, thetas = theta)
probC <- prob_field(polys = gens_shapeC$polys, nrows = n_grid, ncols = n_grid)

#save shape parameters
shapeC <- list("mu" = muC, "kappa" = kappa_C, "sigma" = sigmaC, "C" = C, 
               "theta" = theta)
save(shapeC, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/shapeC.rda")


###################
#shape figures
###################
n_demo <- 3
pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/shapeA.pdf", 
    height = 2, width = 4)
layout(matrix(c(1, 2, 3, 10, 10, 10, 4, 5, 6, 10, 10, 10, 7, 8, 9, 10, 10, 10), 
              byrow = T, ncol = 6))
par(oma=c(0,2,2,5), mar = c(0, 0, 0, 0)) 
for (i in 1:n_demo) {
  plot(gens_shapeA_1$polys[[i]], lwd = .35, xlim = c(.1, .9), ylim = c(.1, .9))
}
for (i in 1:n_demo) {
  plot(gens_shapeA_2$polys[[i]], lwd = .35, xlim = c(.1, .9), ylim = c(.1, .9))
}
for (i in 1:n_demo) {
  plot(gens_shapeA_4$polys[[i]], lwd = .35, xlim = c(.1, .9), ylim = c(.1, .9))
}
mtext("Shape A", side = 3, outer = TRUE, font = 1, line = 0)
mtext(expression(paste(kappa, " = 1 ")), side = 2, at = .85, outer = TRUE, 
      cex = .75)
mtext(expression(paste(kappa, " = 2 ")), side = 2, at = .5, outer = TRUE,
      cex = .75)
mtext(expression(paste(kappa, " = 5 ")), side = 2, at = .15, outer = TRUE,
      cex = .75)
par(mar = rep(1, 4))
image(probA_1,  zlim = c(0, 1),xlim = c(0, 1), ylim = c(0, 1),
      xaxt = "n", yaxt = "n", col = viridis(10))
par(oma=c( 0,0,2,1))
image.plot(legend.only=TRUE, zlim = c(0,1), col = viridis(10)) 
dev.off()

pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/shapeB.pdf", 
    height = 1, width = 3)
n_demo <- 3
par(oma=c(.5,0,2,2), mar = c(0, 0, 0, 1), mfrow = c(1, 4)) 
for (i in 1:n_demo) {
  plot(gens_shapeB$polys[[i]],  lwd = 1, xlim = c(.1, .9), ylim = c(.1, .9))
}
image.plot(probB, zlim = c(0, 1), xlim = c(0, 1), ylim = c(0, 1),
      xaxt = "n", yaxt = "n", col = viridis(10))
mtext("Shape B", side = 3, outer = TRUE, font = 1, line = 0.5)
dev.off()


pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/shapeC.pdf", 
    height = 1, width = 3)
n_demo <- 3
par(oma=c(.5,0,2,2), mar = c(0, 0, 0, 1), mfrow = c(1, 4)) 
for (i in 1:n_demo) {
  plot(gens_shapeC$polys[[i]],  lwd = 1, xlim = c(.1, .9), ylim = c(.1, .9))
}
image.plot(probC, zlim = c(0, 1), xlim = c(0, 1), ylim = c(0, 1),
           xaxt = "n", yaxt = "n", col = viridis(10))
mtext("Shape C", side = 3, outer = TRUE, font = 1, line = 0.5)
dev.off()

