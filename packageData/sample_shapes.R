rm(list = ls())
library("viridis")
library("fields")
library("ContouR")

##############################
#constants across shapes
##############################
n_grid <- 300
n_gen <- 100
p <- 30
Cx <- .5; Cy <- .5
theta_space <- 2*pi/p
theta <- seq(theta_space/2, 2*pi, by = theta_space)

###################
#stop sign (convex)
###################
#write down coordinates of roughly-desired shape (a stop sign)
template_stop <- rbind(c(.6, .5), c(.6, .55), c(.575, .575), c(.55, .6), 
                      c(.5, .6), c(.45, .6), c(.425, .575), c(.4, .55), 
                      c(.4, .5), c(.4, .45), c(.425, .425), c(.45, .4), 
                      c(.5, .4), c(.55, .4), c(.575, .425), c(.6, .45))
#re-scale to take up whole space
template_stop <- rescale(coords = template_stop, eps = .1)

#make poly and find means
stop_poly <- make_poly(template_stop, "stop")
mu_stop <- paral_lengths(p, stop_poly, c(Cx, Cy))

#other parameters
kappa_stop <- 1
sigma_stop <- c(rep(.029, 10), rep(.024, 10), rep(.029, 5), rep(.022, 5))

#generate probability distribution
gens_stop <- gen_conts(n_sim = n_gen, mu_stop, kappa = kappa_stop,
                       sigma_stop, Cx = Cx, Cy = Cy, theta)
stop_prob <- prob_field(polys = gens_stop$polys, nrows = n_grid, ncols = n_grid)

#save shape parameters
stop_sign <- list("mu" = mu_stop, "kappa" = kappa_stop, "sigma" = sigma_stop, 
                  "Cx" = Cx, "Cy" = Cy, "theta" = theta)
save(stop_sign, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/stop_sign.rda")


#####################
#bow tie (non-convex)
#####################
#write down coordinates of roughly-desired shape (a bow tie)
template_tie <- rbind(c(.6, .5), c(.6, .6), c(.5, .55),
                      c(.4, .6), c(.4, .5), c(.4, .4),
                      c(.5, .45), c(.6, .4))

#re-scale to take up whole space
template_tie <- rescale(coords = template_tie, eps = .1)

#make poly and find means
tie_poly <- make_poly(template_tie, "tie")
mu_tie <- paral_lengths(p, tie_poly, c(Cx, Cy))

#other parameters
kappa_tie <- 2
sigma_tie <- c(seq(.025, .033, length = 8), seq(.025, .033, length = 7),
                seq(.033, .025, length = 7), seq(.033, .025, length = 8))

#generate probability distribution
gens_tie <- gen_conts(n_sim = n_gen, mu_tie, kappa = kappa_tie,
                      sigma = sigma_tie, Cx = Cx, Cy = Cy, theta)
tie_prob <- prob_field(polys = gens_tie$polys, nrows = n_grid, ncols = n_grid)


#save shape parameters
tie <- list("mu" = mu_tie, "kappa" = kappa_tie, "sigma" = sigma_tie, 
            "Cx" = Cx, "Cy" = Cy, "theta" = theta)
save(tie, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/tie.rda")


##############################
#star (very non-convex)
##############################
#write down coordinates of roughly-desired shape (a stop sign)
template_tree <- rbind(c(.52, .52), c(.64, .58), c(.51, .54),
                       c(.5, .6), c(.49, .54), c(.36, .58),
                       c(.48, .53),c(.36, .52), c(.48, .52), 
                       c(.36, .48), c(.48, .5),
                       c(.48, .43), c(.52, .43), c(.52, .5),
                       c(.64, .48), c(.52, .51), c(.64, .52),
                       c(.52, .52))

#re-scale to take up whole space
template_tree <- rescale(coords = template_tree, eps = .1)


#make poly and find means
tree_poly <- make_poly(template_tree, "tree")
mu_tree <- paral_lengths(p, tree_poly, c(Cx, Cy))


#other parameters
kappa_tree <- 1
sigma_tree <- c(rep(.033, 5), rep(.024, 5), rep(.033, 10), rep(.02, 5), 
                rep(.025, 5))
                

#generate probability distribution
gens_tree <- gen_conts(n_sim = n_gen, mu_tree, kappa = kappa_tree,
                       sigma = sigma_tree,  Cx = Cx, Cy = Cy, theta)
tree_prob <- prob_field(polys = gens_tree$polys, nrows = n_grid, ncols = n_grid)

#save shape parameters
tree <- list("mu" = mu_tree, "kappa" = kappa_tree, "sigma" = sigma_tree, 
            "Cx" = Cx, "Cy" = Cy, "theta" = theta)
save(tree, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/tree.rda")


###################
#shape figures
###################
pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/sample_shapes_prob.pdf", 
    height = 4, width = 8.5)
set.panel()
par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel(1,3) 
image(stop_prob, zlim = c(0, 1),
           xaxt = "n", yaxt = "n", col = viridis(10),
           main = "Stop Sign")
image(tie_prob,  zlim = c(0, 1),
           xaxt = "n", yaxt = "n", col = viridis(10),
           main = "Bow Tie")
image(tree_prob, zlim = c(0, 1),
           xaxt = "n", yaxt = "n", col = viridis(10),
           main = "Tree")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim = c(0,1), col = viridis(10)) 
dev.off()

n_demo <- 15
pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/sample_shapes_contours.pdf", height = 4, width = 8.5)
set.panel(1,3) 
plot(gens_stop$polys[[1]], lwd = .25, main = "Stop Sign", xlim = c(0, 1), 
     ylim = c(0, 1))
for (i in 2:n_demo) {
  plot(gens_stop$polys[[i]], add = T, lwd = .35)
}
plot(gens_tie$polys[[1]], lwd = .25, main = "Bow Tie", xlim = c(0, 1), 
     ylim = c(0, 1))
for (i in 2:n_demo) {
  plot(gens_tie$polys[[i]], add = T, lwd = .35)
}
plot(gens_tree$polys[[1]],  lwd = .25, main = "Tree", xlim = c(0, 1), 
     ylim = c(0, 1))
for (i in 2:n_demo) {
  plot(gens_tree$polys[[i]], add = T, lwd = .35)
}
dev.off()




