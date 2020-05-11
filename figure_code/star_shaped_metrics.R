rm(list = ls())
library("ContouR")
library("maptools")
library("xtable")
p <- 200
theta_space <- 2*pi/p
thetas <- seq(theta_space/2, 2*pi, theta_space)
box <- bbox()
bd <- box@polygons[[1]]@Polygons[[1]]@coords

###################################
#Compute metrics for the 4 shapes
#####################################
cp_poly <- list()
cp_poly[[1]] <- make_poly(convex_poly, "cp")
cp_C_under <- cp_C_over <- matrix(data = c(.5, .5), ncol = 2)
cp <- assess_star(cont = cp_poly[[1]], C_under = cp_C_under, C_over = cp_C_over,
                  thetas = thetas)
cm_poly <- list()
cm_poly[[1]] <- make_poly(convex_mod, "cm")
cm_C_under <- best_C(bd = bd, conts = cm_poly, thetas = thetas, space = .04)
cm_C_over <- best_C(bd = bd, conts = cm_poly, thetas = thetas,  space = .04,
                    under = FALSE)
cm <- assess_star(cont = cm_poly[[1]], C_under = cm_C_under, C_over = cm_C_over,
                  thetas = thetas)

circl_poly <- list()
circl_poly[[1]] <- make_poly(circle_mod, "circle")
circl_C_under <- best_C(bd = bd, conts = circl_poly, thetas = thetas, 
                        space = .05)
circl_C_over <- best_C(bd = bd, conts = circl_poly, thetas = thetas, 
                       space  = .05, under = FALSE)
circl <- assess_star(cont = circl_poly[[1]], C_under = circl_C_under, 
                     C_over = circl_C_over, thetas = thetas)

curl_poly <- list()
curl_poly[[1]] <- make_poly(curl, "circle")
curl_C_under <- best_C(bd = bd, conts = curl_poly, thetas = thetas, 
                       space = .05)
curl_C_over <- best_C(bd = bd, conts = curl_poly, thetas = thetas, 
                      space = .05, under = FALSE)
curl <- assess_star(cont = curl_poly[[1]], C_under = curl_C_under, 
                    C_over = curl_C_over, thetas = thetas)

#########################
#make figure
#########################
pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/starShapedMetric.pdf")
par(mfrow = c(4, 3), oma = c(0, .5, 2.8, 0), mar = rep(0, 4))
#ex 1
plot(cp_poly[[1]])
plot(cp_poly[[1]])
plot(cp$u_approx, add = T, border = 'red')
points(cp_C_under, pch = 3, col = 'red', cex = 2)
plot(cp_poly[[1]])
plot(cp$o_approx, add = T, border = 'blue')
points(cp_C_over, pch = 3, col = 'blue', cex = 2)

#example 2
plot(cm_poly[[1]])
plot(cm_poly[[1]])
plot(cm$u_error, add = T, col = 'pink')
plot(cm$u_approx, add = T, border = 'red')
points(cm_C_under, pch = 3, col = 'red', cex = 2)
plot(cm_poly[[1]])
plot(cm$o_error, add = T, col = 'lightblue')
plot(cm$o_approx, add = T, border = 'blue')
points(cm_C_over, pch = 3, col = 'blue', cex = 2)

#example 3
plot(circl_poly[[1]])
plot(circl_poly[[1]])
plot(circl$u_error, add = T, col = 'pink')
plot(circl$u_approx, add = T, border = 'red')
points(circl_C_under, pch = 3, col = 'red', cex = 2)
plot(circl_poly[[1]])
plot(circl$o_error, add = T, col = 'lightblue')
plot(circl$o_approx, add = T, border = 'blue')
points(circl_C_over, pch = 3, col = 'blue', cex = 2)

#example 4
plot(curl_poly[[1]])
plot(curl_poly[[1]])
plot(curl$u_error, add = T, col = 'pink')
points(curl_C_under, pch = 3, col = 'red')
plot(curl$u_approx, add = T, border = 'red', cex = 2)
plot(curl_poly[[1]])
plot(curl$o_error, add = T, col = 'lightblue')
plot(curl$o_approx, add = T, border = 'blue')
points(curl_C_over, pch = 3, col = 'blue', cex = 2)

#labels
mtext(text = "Contour", side = 3, line = 0, at = .18, outer = TRUE)
mtext(text = "Under-estimated.", side = 3, line = 0, at = .5, 
      outer = TRUE)
mtext(text = "Over-estimated", side = 3, line = 0, at = .82,
      outer = TRUE)
mtext(text = "Contour 1", side = 2, line = -2, at = .89, outer = TRUE)
mtext(text = "Contour 2", side = 2, line = -2, at = .63, outer = TRUE)
mtext(text = "Contour 3", side = 2, line = -2, at = .38, outer = TRUE)
mtext(text = "Contour 4", side = 2, line = -2, at = .1, outer = TRUE)
dev.off()

######################
#table of areas
#####################
prop_areas <- matrix(nrow = 4, ncol = 2, byrow = TRUE,
                      data = c(cp$u_area/cp$tot_area, cp$o_area/cp$tot_area,
                               cm$u_area/cm$tot_area, cm$o_area/cm$tot_area,
                               circl$u_area/circl$tot_area, circl$o_area/circl$tot_area,
                               curl$u_area/curl$tot_area, curl$o_area/circl$tot_area))
per_areas <- 100*prop_areas
xtable(per_areas, digits = 2)
