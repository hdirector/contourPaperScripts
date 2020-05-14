rm(list = ls())
library("ContouR")
library("xtable")
p <- 200
theta_space <- 2*pi/p
thetas <- seq(theta_space/2, 2*pi, theta_space)
box <- bbox()
bd <- box@polygons[[1]]@Polygons[[1]]@coords

###################################
#Compute metrics for the 4 shapes
#####################################
#shape 1
cp_poly <- list()
cp_poly[[1]] <- make_poly(convex_poly, "cp")
cp_C_under <- cp_C_over <- matrix(data = c(.5, .5), ncol = 2)
cp_over <- error_areas(conts = cp_poly, C = cp_C_under, theta = thetas,
                       under  = FALSE, ret_all = TRUE)
cp_under <- error_areas(conts = cp_poly, C = cp_C_over, theta = thetas,
                        under  = TRUE, ret_all = TRUE)

#shape 2
cm_poly <- list()
cm_poly[[1]] <- make_poly(convex_mod, "cm")
cm_C_under <- best_C(bd = bd, conts = cm_poly, thetas = thetas, space = .04,
                     parallel = TRUE)
cm_C_over <- best_C(bd = bd, conts = cm_poly, thetas = thetas,  space = .04,
                    under = FALSE, parallel = TRUE)
cm_over <- error_areas(conts = cm_poly, C = cm_C_over, theta = thetas,
                      under  = FALSE, ret_all = TRUE)
cm_under <- error_areas(conts = cm_poly, C = cm_C_under, theta = thetas,
                       under  = TRUE, ret_all = TRUE)

#shape 3
circl_poly <- list()
circl_poly[[1]] <- make_poly(circle_mod, "circle")
circl_C_under <- best_C(bd = bd, conts = circl_poly, thetas = thetas, 
                        space = .05, parallel = TRUE)
circl_C_over <- best_C(bd = bd, conts = circl_poly, thetas = thetas, 
                       space  = .05, under = FALSE, parallel = TRUE)
circl_over <- error_areas(conts = circl_poly, C = circl_C_over, theta = thetas,
                          under  = FALSE, ret_all = TRUE)
circl_under <- error_areas(conts = circl_poly, C = circl_C_under, theta = thetas,
                           under  = TRUE, ret_all = TRUE)

#shape 4
curl_poly <- list()
curl_poly[[1]] <- make_poly(curl, "circle")
curl_C_under <- best_C(bd = bd, conts = curl_poly, thetas = thetas, 
                       space = .05, parallel = TRUE)
curl_C_over <- best_C(bd = bd, conts = curl_poly, thetas = thetas, 
                      space = .05, under = FALSE, parallel = TRUE)
curl_over <- error_areas(conts = curl_poly, C = curl_C_over, theta = thetas,
                         under  = FALSE, ret_all = TRUE)
curl_under <- error_areas(conts = curl_poly, C = curl_C_under, theta = thetas,
                          under  = TRUE, ret_all = TRUE)

#########################
#make figure
#########################
pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/starShapedMetric.pdf")
par(mfrow = c(4, 3), oma = c(0, .5, 2.8, 0), mar = rep(0, 4))
#shape 1
plot(cp_poly[[1]])
plot(cp_poly[[1]])
plot(cp_under$approxs[[1]], add = T, border = 'red')
points(cp_C_under, pch = 3, col = 'red', cex = 2)
plot(cp_poly[[1]])
plot(cp_over$approxs[[1]], add = T, border = 'blue')
points(cp_C_over, pch = 3, col = 'blue', cex = 2)

#example 2
plot(cm_poly[[1]])
plot(cm_poly[[1]])
plot(cm_under$diffs[[1]], add = T, col = 'pink')
plot(cm_under$approxs[[1]], add = T, border = 'red')
points(cm_C_under, pch = 3, col = 'red', cex = 2)
plot(cm_poly[[1]])
plot(cm_over$diffs[[1]], add = T, col = 'lightblue')
plot(cm_over$approxs[[1]], add = T, border = 'blue')
points(cm_C_over, pch = 3, col = 'blue', cex = 2)


#example 3
plot(circl_poly[[1]])
plot(circl_poly[[1]])
plot(circl_under$diffs[[1]], add = T, col = 'pink')
plot(circl_under$approxs[[1]], add = T, border = 'red')
points(circl_C_under, pch = 3, col = 'red', cex = 2)
plot(circl_poly[[1]])
plot(circl_over$diffs[[1]], add = T, col = 'lightblue')
plot(circl_over$approxs[[1]], add = T, border = 'blue')
points(circl_C_over, pch = 3, col = 'blue', cex = 2)

#example 4
plot(curl_poly[[1]])
plot(curl_poly[[1]])
plot(curl_under$diffs[[1]], add = T, col = 'pink')
points(curl_C_under, pch = 3, col = 'red', cex = 2)
plot(curl_under$approxs[[1]], add = T, border = 'red', cex = 2)
plot(curl_poly[[1]])
plot(curl_over$diffs[[1]], add = T, col = 'lightblue')
plot(curl_over$approxs[[1]], add = T, border = 'blue')
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
                      data = c(cp_under$areas[[1]]/gArea(cp_poly[[1]]), 
                               cp_over$areas[[1]]/gArea(cp_poly[[1]]),
                               cm_under$areas[[1]]/gArea(cm_poly[[1]]), 
                               cm_over$areas[[1]]/gArea(cm_poly[[1]]),
                               circl_under$areas[[1]]/gArea(circl_poly[[1]]), 
                               circl_over$areas[[1]]/gArea(circl_poly[[1]]),
                               curl_under$areas[[1]]/gArea(curl_poly[[1]]),
                               curl_over$areas[[1]]/gArea(curl_poly[[1]])))
per_areas <- 100*prop_areas
xtable(per_areas, digits = 2)
