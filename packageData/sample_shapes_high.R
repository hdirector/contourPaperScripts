rm(list = ls())
library("viridis")
library("fields")
library("ContouR")

##############################
#constants across shapes
##############################
n_grid <- 300
n_gen <- 100
p <- 100
Cx <- .5; Cy <- .5
theta_space <- 2*pi/p
theta <- seq(theta_space/2, 2*pi, by = theta_space)

###################
#stop sign (convex)
###################
#write down coordinates of roughly-desired shape (a stop sign)
theta_gen <- seq(2*pi/20, 2*pi, by = 2*pi/20)
template_stop_hd <- cbind(.5 + .1*cos(theta_gen), .5 +.1*sin(theta_gen))
stop_hd_poly <- make_poly(template_stop_hd, "stop_hd")
mu_stop_hd <- paral_lengths(p, stop_hd_poly, c(Cx, Cy))

#other parameters
kappa_stop_hd <- 1
sigma_stop_hd <- c(rep(.009, 5), rep(.005, 10), rep(.012, 10), rep(.007, 5),
                   rev(c(rep(.009, 5), rep(.005, 10), rep(.012, 10), rep(.007, 5))),
                   rep(.009, 5), rep(.005, 10), rep(.012, 10), rep(.007, 5),
                   rep(.012, 10))

#generate probability distribution
gens_stop_hd <- gen_conts(n_sim = n_gen, mu_stop_hd, kappa = kappa_stop_hd,
                       sigma_stop_hd, Cx = Cx, Cy = Cy, theta)
stop_hd_bound <- sapply(gens_stop_hd$polys, function(x){conv_to_grid(x, n_grid, n_grid)})
stop_hd_prob <- matrix(apply(stop_hd_bound, 1, mean), nrow = n_grid, ncol = n_grid)

#save shape parameters
stop_sign_hd <- list("mu" = mu_stop_hd, "kappa" = kappa_stop_hd, "sigma" = sigma_stop_hd, 
                  "Cx" = Cx, "Cy" = Cy, "theta" = theta)
save(stop_sign_hd, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/stop_sign_hd.rda")


#####################
#bow tie (non-convex)
#####################
#write down coordinates of roughly-desired shape (a bow tie)
template_tie_hd <- rbind(c(.6, .5), c(.6, .54), c(.58, .55), c(.6, .56), 
                      c(.6, .6), c(.56, .58), c(.55, .565), c(.54, .57), c(.5, .55),
                      c(.46, .57), c(.45, .565), c(.44, .58), 
                      c(.4, .6), c(.4, .56), c(.42, .55),
                      c(.4, .54), 
                      c(.4, .5), c(.4, .46), c(.42, .45), c(.4, .44),  
                      c(.4, .4), c(.44, .42), c(.45, .435), c(.46, .43), 
                      c(.5, .45), c(.54, .43), c(.55, .435), c(.56, .42), 
                      c(.6, .4), c(.6, .44), c(.58, .45), c(.6, .46))
tie_hd_poly <- make_poly(template_tie_hd, "tie_hd")

#convert coordinates to parameters used in demos
mu_tie_hd <- paral_lengths(p, tie_hd_poly, c(Cx, Cy))

#other parameters
kappa_tie_hd <- 2
sigma_tie_hd <- c(rep(c(seq(.007, .014, length = 8), seq(.014, .007, length = 7),
                  seq(.007, .014, length = 7), seq(.014, .007, length = 8)), 3),
               rep(.012, 10))

#generate probability distribution
gens_tie_hd <- gen_conts(n_sim = n_gen, mu_tie_hd, kappa = kappa_tie_hd,
                      sigma = sigma_tie_hd, Cx = Cx, Cy = Cy, theta)
tie_hd_bound <- sapply(gens_tie_hd$polys, function(x){conv_to_grid(x, n_grid, n_grid)})
tie_hd_prob <- matrix(apply(tie_hd_bound, 1, mean), nrow = n_grid, ncol = n_grid)


#save shape parameters
tie_hd <- list("mu" = mu_tie_hd, "kappa" = kappa_tie_hd, "sigma" = sigma_tie_hd, 
            "Cx" = Cx, "Cy" = Cy, "theta" = theta)
save(tie_hd, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/tie_hd.rda")


##############################
#star (very non-convex)
##############################
#write down coordinates of roughly-desired shape (a stop sign)
template_tree_hd <- rbind(c(.52, .52), c(.64, .58), c(.51, .54),
                       c(.5, .6), c(.49, .54), c(.36, .58),
                       c(.48, .53),c(.36, .52), c(.48, .52), 
                       c(.36, .48), c(.48, .5),
                       c(.48, .43), c(.52, .43), c(.52, .5),
                       c(.64, .48), c(.52, .51), c(.64, .52),
                       c(.52, .52))


tree_hd_poly <- make_poly(template_tree_hd, "tree_hd")

#convert coordinates to parameters used in demos
mu_tree_hd <- paral_lengths(p, tree_hd_poly, c(Cx, Cy))


#other parameters
kappa_tree_hd <- 1
sigma_tree_hd <- c(rep(c(rep(.011, 5), rep(.008, 5), rep(.012, 10), rep(.008, 5), 
                   rep(.009, 5)), 3), rep(.012, 10))


#generate probability distribution
gens_tree_hd <- gen_conts(n_sim = n_gen, mu_tree_hd, kappa = kappa_tree_hd,
                       sigma = sigma_tree_hd,  Cx = Cx, Cy = Cy, theta)
tree_hd_bound <- sapply(gens_tree_hd$polys, function(x){conv_to_grid(x, n_grid, n_grid)})
tree_hd_prob <- matrix(apply(tree_hd_bound, 1, mean), nrow = n_grid, ncol = n_grid)

#save shape parameters
tree_hd <- list("mu" = mu_tree_hd, "kappa" = kappa_tree_hd, "sigma" = sigma_tree_hd, 
             "Cx" = Cx, "Cy" = Cy, "theta" = theta)
#save(tree, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/tree_hd.rda")


###################
#shape figures
###################
#pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/sample_shapes_hd_prob.pdf", 
#    height = 4, width = 8.5)
set.panel()
par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel(1,3) 
image(stop_hd_prob, xlim = c(.3, .7), ylim = c(.3, .7), zlim = c(0, 1),
      xaxt = "n", yaxt = "n", col = viridis(10),
      main = "Stop Sign")
image(tie_hd_prob, xlim = c(.3, .7), ylim = c(.3, .7), zlim = c(0, 1),
      xaxt = "n", yaxt = "n", col = viridis(10),
      main = "Bow Tie")
image(tree_hd_prob, xlim = c(.3, .7), ylim = c(.3, .7), zlim = c(0, 1),
      xaxt = "n", yaxt = "n", col = viridis(10),
      main = "Tree")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim = c(0,1), col = viridis(10)) 
#dev.off()

n_demo <- 15
#pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/sample_shapes_hd_contours.pdf", height = 4, width = 8.5)
set.panel(1,3) 
plot(gens_stop_hd$polys[[1]], xlim = c(.3, .7), ylim = c(.3, .7),
     lwd = .25, main = "Stop Sign")
for (i in 2:n_demo) {
  plot(gens_stop_hd$polys[[i]], add = T, lwd = .35)
}
plot(gens_tie_hd$polys[[1]], xlim = c(.3, .7), ylim = c(.3, .7),
     lwd = .25, main = "Bow Tie")
for (i in 2:n_demo) {
  plot(gens_tie_hd$polys[[i]], add = T, lwd = .35)
}
plot(gens_tree_hd$polys[[1]], xlim = c(.3, .7), ylim = c(.3, .7),
     lwd = .25, main = "Tree")
for (i in 2:n_demo) {
  plot(gens_tree_hd$polys[[i]], add = T, lwd = .35)
}
#dev.off()




