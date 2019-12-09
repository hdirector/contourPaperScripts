rm(list = ls())
#TO DO: ADD FINDING THE OPTIMAL C EACH EXAMPLE
library("ContouR")
library("maptools")
library("xtable")
p <- 500
theta <- seq(0, 2*pi, length = p)

###################################
#Compute metrics for the 4 shapes
#####################################
cp_poly <- make_poly(convex_poly, "cp")
cp <- assess_star(cont = cp_poly, C = c(.5, .5),theta = theta)
cm_poly <- make_poly(convex_mod, "cm")
cm <- assess_star(cont = cm_poly, C = c(.5, .5),theta = theta)
circl_poly <- make_poly(circle_mod, "circle")
circl <- assess_star(cont = circl_poly, C = c(.5, .5),theta = theta)
curl_poly <- make_poly(curl, "circle")
curl <- assess_star(cont = curl_poly, C = c(.5, .5),theta = theta)


###########################
#make figure
#########################
pdf("/users/hdirector/desktop/starShapedDemo.pdf")
par(mfrow = c(4, 3), oma = rep(.2, 4), mar = c(.2, .2, 4, .2))
#ex 1
plot(cp_poly)
plot(cp_poly)
plot(cp$u_error, add = T, col = 'pink')
plot(cp$u_approx, add = T, border = 'red')
plot(cp_poly)
plot(cp$o_error, add = T, col = 'lightblue')
plot(cp$o_approx, add = T, border = 'blue')

#example 2
plot(cm_poly)
plot(cm_poly)
plot(cm$u_error, add = T, col = 'pink')
plot(cm$u_approx, add = T, border = 'red')
plot(cm_poly)
plot(cm$o_error, add = T, col = 'lightblue')
plot(cm$o_approx, add = T, border = 'blue')

#example 3
plot(circl_poly)
plot(circl_poly)
plot(circl$u_error, add = T, col = 'pink')
plot(circl$u_approx, add = T, border = 'red')
plot(circl_poly)
plot(circl$o_error, add = T, col = 'lightblue')
plot(circl$o_approx, add = T, border = 'blue')

#example 4
plot(curl_poly)
plot(curl_poly)
plot(curl$u_error, add = T, col = 'pink')
plot(curl$u_approx, add = T, border = 'red')
plot(curl_poly)
plot(curl$o_error, add = T, col = 'lightblue')
plot(curl$o_approx, add = T, border = 'blue')

mtext(text = "Contour", side = 3, line = -3, at = .18, outer = TRUE)
mtext(text = "Under Approx \n Error Area", side = 3, line = -3, at = .5, 
      outer = TRUE)
mtext(text = "Over Approx. \n Error Area", side = 3, line = -3, at = .82,
      outer = TRUE)
mtext(text = "Contour 1", side = 2, line = -3, at = .85, outer = TRUE)
mtext(text = "Contour 2", side = 2, line = -3, at = .58, outer = TRUE)
mtext(text = "Contour 3", side = 2, line = -3, at = .35, outer = TRUE)
mtext(text = "Contour 4", side = 2, line = -3, at = .1, outer = TRUE)
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
xtable(per_areas, digits = 3)