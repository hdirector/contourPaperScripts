rm(list = ls())
#TO DO: ADD FINDING THE OPTIMAL C EACH EXAMPLE
library("ContouR")
library("maptools")
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



####################
#convex mod
####################
cm <- make_poly(convex_mod, "cm_poly")
C_cm <- c(.5, .5)
l_cm <- make_l(C_cm, theta, r = 5)
u_cm <- make_poly(pts_on_l(l = l_cm, cont = cm, under = TRUE), "u_cm")
o_cm <- make_poly(pts_on_l(l = l_cm, cont = cm, under = FALSE), "o_cm")

u_area <- diff_reg(cm, u_cm)
o_area <- diff_reg(u_cm, cm) 
gArea(u_area)
gArea(o_area)

#plot to demonstrate
par(mfrow = c(1, 3))
plot(cm, main = "Contour")
plot(cm, main = "Under Approximation")
plot(u_area, add = TRUE, col = 'pink')
plot(u_cm, add = TRUE, border = 'red')
plot(cm, main = "Over Approximation")
plot(o_cm, add = TRUE, border = 'blue')
plot(o_area, add = TRUE, col = 'lightblue')
plot(o_cm, add = TRUE, border = 'blue')

####################
#circle mod
####################

####################
#curl
####################
assess_star(cur)
curl <- make_poly(curl, "curl_poly")
C_curl <- c(.5, .5)
l_curl <- make_l(C_curl, theta, r = 5)
u_curl <- make_poly(pts_on_l(l = l_curl, cont = curl, under = TRUE), "u_curl")
o_curl <- make_poly(pts_on_l(l = l_curl, cont = curl, under = FALSE), "o_curl")

u_area <- diff_reg(curl, u_curl)
o_area <- diff_reg(u_curl, curl) 
gArea(u_area)
gArea(o_area)

#plot to demonstrate
par(mfrow = c(1, 3))
plot(curl, main = "Contour")
plot(curl, main = "Under Approximation")
plot(u_area, add = TRUE, col = 'pink')
plot(u_curl, add = TRUE, border = 'red')
plot(curl, main = "Over Approximation")
plot(o_curl, add = TRUE, border = 'blue')
plot(o_area, add = TRUE, col = 'lightblue')
plot(o_curl, add = TRUE, border = 'blue')


