library("ContouR")
##############################
#constants across shapes
##############################
n_grid <- 150
n_gen <- 100
p <- 50
Cx <- .5; Cy <- .5
theta_space <- 2*pi/p
theta <- seq(theta_space/2, 2*pi, by = theta_space)

##########################
# parameter settings
##########################
circ <- list(C = c(0.5, 0.5), a = 0, b = .3, d = 0, 
             k = .03, l = .01, m = 1,
             kappa = .3)
fig8 <- list(C = c(0.5, 0.5), a = .15, b = .1, d = 2, 
             k = .01, l = .025, m = 2, 
             kappa = .2)
flower <- list(C = c(0.5, 0.5), a = .16, b = .08, d = 6, 
               k = .03, l = .03, m = 0, 
               kappa = 1.2)

save(circ, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/circ.rda")
save(fig8, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/fig8.rda")
save(flower, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/flower.rda")

##############################
#sampling settings
##############################
n_grid <- 150
n_gen <- 100
theta_space <- 2*pi/p
theta <- seq(theta_space/2, 2*pi, by = theta_space)

#Circle
pts_circ <- discrete_pts(circ$C, circ$a, circ$b, circ$d, theta)
mu_circ <- discrete_mu(circ$a, circ$b, circ$d, theta)
sigma_circ <- discrete_sigma(circ$k, circ$l, circ$m, theta)

#figure eight
pts_fig8 <- discrete_pts(fig8$C, fig8$a, fig8$b, fig8$d, theta)
mu_fig8 <- discrete_mu(fig8$a, fig8$b, fig8$d, theta)
sigma_fig8 <- discrete_sigma(fig8$k, fig8$l, fig8$m, theta)

#flower
pts_flower <- discrete_pts(flower$C, flower$a, flower$b, flower$d,theta)
mu_flower <- discrete_mu(flower$a, flower$b, flower$d, theta)
sigma_flower <- discrete_sigma(flower$k, flower$k, flower$m, theta)

#plot means
par(mfrow = c(1, 3))
plot(pts_circ, type= "l", xlim = c(0, 1), ylim = c(0, 1))
plot(pts_fig8, type= "l", xlim = c(0, 1), ylim = c(0, 1))
plot(pts_flower, type= "l", xlim = c(0, 1), ylim = c(0, 1))


###################################
#generate probability distribution
###################################
gens_circ <- gen_conts(n_sim = n_gen, mu = mu_circ, kappa = circ$kappa,
                       sigma = sigma_circ, C = circ$C, thetas = theta)
circ_prob <- prob_field(polys = gens_circ$polys, nrows = n_grid, ncols = n_grid)

gens_fig8 <- gen_conts(n_sim = n_gen, mu = mu_fig8, kappa = fig8$kappa,
                       sigma = sigma_fig8, C = fig8$C, thetas = theta)
fig8_prob <- prob_field(polys = gens_fig8$polys, nrows = n_grid, ncols = n_grid)

gens_flower <- gen_conts(n_sim = n_gen, mu = mu_flower, kappa = flower$kappa,
                         sigma = sigma_flower, C = flower$C,  thetas = theta)
flower_prob <- prob_field(polys = gens_flower$polys, nrows = n_grid, ncols = n_grid)


###################
#shape figures
###################
#pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/sample_shapes_prob.pdf", 
#    height = 4, width = 8.5)
library("viridis")
library("fields")
set.panel()
par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel(1,3) 
image(circ_prob, zlim = c(0, 1),
      xaxt = "n", yaxt = "n", col = viridis(10),
      main = "Circle")
image(fig8_prob,  zlim = c(0, 1),
      xaxt = "n", yaxt = "n", col = viridis(10),
      main = "Figure 8")
image(flower_prob, zlim = c(0, 1),
      xaxt = "n", yaxt = "n", col = viridis(10),
      main = "Flower")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim = c(0,1), col = viridis(10)) 
#dev.off()

n_demo <- 15
#pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/sample_shapes_contours.pdf", height = 4, width = 8.5)
set.panel(1,3) 
plot(gens_circ$polys[[1]], lwd = .25, main = "Circle", xlim = c(0, 1), 
     ylim = c(0, 1))
for (i in 2:n_demo) {
  plot(gens_circ$polys[[i]], add = T, lwd = .35)
}
plot(gens_fig8$polys[[1]], lwd = .25, main = "Figure 8", xlim = c(0, 1), 
     ylim = c(0, 1))
for (i in 2:n_demo) {
  plot(gens_fig8$polys[[i]], add = T, lwd = .35)
}
plot(gens_flower$polys[[1]],  lwd = .25, main = "Flower", xlim = c(0, 1), 
     ylim = c(0, 1))
for (i in 2:n_demo) {
  plot(gens_flower$polys[[i]], add = T, lwd = .35)
}
#dev.off()


