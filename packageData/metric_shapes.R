#make points for four shapes to demonstrate the value of metrics
rm(list = ls())
library("ContouR")
p <- 200
theta_space <- 2*pi/p
thetas <- seq(theta_space/2, 2*pi, by = theta_space)
C <- c(.5, .5)

##########################
#Star-shaped (stop sign)
##########################
cp_pts_orig <- rbind(c(.6, .5), c(.6, .55), c(.575, .575), c(.55, .6), 
                  c(.5, .6), c(.45, .6), c(.425, .575), c(.4, .55), 
                  c(.4, .5), c(.4, .45), c(.425, .425), c(.45, .4), 
                  c(.5, .4), c(.55, .4), c(.575, .425), c(.6, .45))
cp_pts_orig <- rbind(cp_pts_orig, cp_pts_orig[1,])
#re-scale 
cp_pts <- rescale(coords = list(cp_pts_orig), eps = .1, 
                  box_size = 1)$coords_scale[[1]]
cp_mu <- paral_lengths(p = p, poly = make_poly(convex_poly, "convex_poly"), 
                     C = C)
convex_poly <- cbind(C[1] + cp_mu*cos(thetas), C[2] + cp_mu*sin(thetas))


##########################
#Star-shaped (stop sign)
##########################
cm_pts_orig <- rbind(c(.6, .5), c(.6, .55), c(.575, .575), c(.62, .62),
                     c(.58, .61), c(.58, .59),c(.55, .6), 
                     c(.5, .6), c(.45, .6), c(.425, .575), c(.4, .55), 
                     c(.4, .5), c(.4, .45), c(.425, .425), c(.45, .4), 
                     c(.5, .4), c(.55, .4), c(.575, .425), c(.6, .45))
cm_pts_orig <- rbind(cm_pts_orig, cm_pts_orig[1,])
#re-scale 
convex_mod_poly <- rescale(coords = list(cm_pts_orig), eps = .1,
                           box_size = 1)$coords_scale[[1]]


#####################
#modified circle
######################
r1 <- .3; r2 <- .55
theta1 <- seq(pi/10, pi/4, length = 3)
theta2 <- seq(pi/4 + pi/10, pi/2, length = 3)
theta3 <- seq(pi/2 + pi/10, 3*pi/4, length = 3)
theta4 <- seq(3*pi/4 + pi/10, pi, length = 3)
theta5 <- seq(pi + pi/10, 5*pi/4, length = 3)
theta6 <- seq(5*pi/4 + pi/10, 3*pi/2, length = 3)
theta7 <- seq(3*pi/2 + pi/10, 7*pi/4, length = 3)
theta8 <- seq(7*pi/4 + pi/10, 2*pi, length = 3)
circle_mod_pts <- rbind(cbind(r1*cos(theta1), r1*sin(theta1)),
                    c(r2*cos(theta1[2]), r2*sin(theta1[2])),
                    cbind(r1*cos(theta2), r1*sin(theta2)),
                    cbind(r1*cos(theta3), r1*sin(theta3)),
                    c(r2*cos(theta3[2]), r2*sin(theta3[2])),
                    cbind(r1*cos(theta4), r1*sin(theta4)),
                    cbind(r1*cos(theta5), r1*sin(theta5)),
                    c(r2*cos(theta5[2]), r2*sin(theta5[2])),
                    cbind(r1*cos(theta6), r1*sin(theta6)),
                    cbind(r1*cos(theta7), r1*sin(theta7)),
                    c(r2*cos(theta7[2]), r2*sin(theta7[2])),
                    cbind(r1*cos(theta8), r1*sin(theta8)))
#re-scale 
circle_mod <- rescale(coords = list(circle_mod_pts), eps = .1, box_size = 1)$coords_scale[[1]]


####################
#curl
#####################
r1 <- .3; r2 <- .5; r3 <- .4
theta1 <- seq(0, 7*pi/4, length = 15)
theta2 <- c(seq(7*pi/4, 2*pi,length = 4), seq(0, pi/2, length = 5))
theta3 <- seq(pi/2, 0, length = 5)
curl_pts <- rbind(cbind(r1*cos(theta1), r1*sin(theta1)),
                  cbind(r2*cos(theta2), r2*sin(theta2)),
                  cbind(r3*cos(theta3), r3*sin(theta3)))
curl_pts <- rbind(curl_pts, curl_pts[1,])
#re-scale 
curl <- rescale(coords = list(curl_pts), eps = .1, box_size = 1)$coords_scale[[1]]

################
#plot shapes
################
par(mfrow = c(1, 4))
plot(convex_poly, type = "l", main = "Convex Polygon", xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n")
plot(convex_mod, type = "l", main = "Modified Stop Sign", xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n")
plot(circle_mod, type = "l", main = "Modified Circle", xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n")
plot(curl, type = "l", main = "Curl", xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n")

##########################
#save and add to package
##########################
save(convex_poly, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/convex_poly.rda")
save(convex_mod, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/convex_mod.rda")
save(circle_mod, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/circle_mod.rda")
save(curl, file = "/Users/hdirector/Dropbox/Contours/ContouR/data/curl.rda")

