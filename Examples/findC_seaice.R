library("IceCast")
library("ContouR")
library("raster")

#read in September sea ice observations
all <- IceCast::read_monthly_BS(1980, 2017, 
                                "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/bootstrapV3_1/",
                                version = 3.1)
sep <- all[,9,,]
n_obs <- dim(sep)[1]

#regions in and around central Arctic
reg <- reg_info$regions[[1]]
reg_info$regions[[1]]@polygons[[1]]@ID <- "centArc"
reg <- spRbind(reg, reg_info$regions[[2]])
reg <- rm_holes(raster::aggregate(reg))

#make a grid of points at the center of each square (
xmn = -3850; xmx = 3750; ymn = -5350; ymx = 5850
x_cent <- seq(xmn + 13.5, xmx, 150)
y_cent <- seq(ymn + 13.5, ymx, 150)
grid <- expand.grid(x_cent, y_cent)

coords <- list()
for (i in 1:n_obs) {
  #find coordinates of largest contigous region in Central Arctic
  temp <- disaggregate(get_region(sep[i,,], dat_type = "bootstrap", level = 15))
  temp <- temp[which.max(gArea(temp, byid = TRUE))]
  temp <- disaggregate(rm_holes(keep_poly(gIntersection(reg, temp))))
  temp <- temp[which.max(gArea(temp, byid = T))]
  coords[[i]] <- temp@polygons[[1]]@Polygons[[1]]@coords
}

#rescale everything 
temp_rescale <- rescale(coords, eps = .1, grid)
coords_scale <- temp_rescale$coords_scale
grid_scale <- temp_rescale$grid_scale

pdf("/users/hdirector/desktop/rescale.pdf", height = 4, width = 8)
par(mfrow = c(1, 2))
plot(coords[[1]], type= "l", xlab = "", ylab = "", main = "1980 Contour")
plot(coords_scale[[1]], type= "l", xlab = "", ylab = "",
     xlim = c(0, 1), ylim = c(0, 1), main = "1980 Contour Rescaled")
dev.off()

#Switch to Spatial objects
grid_pts_temp <- SpatialPoints(grid_scale)
conts <- lapply(coords_scale, function(x){make_poly(x, "contours")})

#Only test points that are in every contour
test <- sapply(conts, function(x){gIntersects(grid_pts_temp, x, byid = TRUE)})
keep <- apply(test, 1, function(x){all(x)})
grid_pts <- SpatialPoints(grid_scale[keep,])

# rbPal <- colorRampPalette(c('red','blue'))
# cols <- rbPal(10)[as.numeric(cut(log(test), breaks = 10))]
# points(grid_pts, col = cols)
# a <- which.min(test)
# points(grid_pts[a,], col = "green")
n_test_pts <- nrow(grid_pts@coords)
p <- 100
r <- 5
theta_dist <- 2*pi/p
theta <- seq(theta_dist/2, 2*pi, by = theta_dist)
C_test <- grid_pts@coords
area_out <- matrix(nrow = n_test_pts, ncol = n_obs)
for (i in 1:n_test_pts) {
  l_pts <- cbind(C_test[i, 1] + r*cos(theta), C_test[i, 2] + r*sin(theta))
  l <- apply(l_pts, 1, function(x){make_line(C_test[i,], x, "l")})
  for (j in 1:n_obs) {
    on_l <- lapply(l, function(x){gIntersection(conts[[j]], x)})
    pts_on_l <- t(sapply(on_l, function(x){x@lines[[1]]@Lines[[1]]@coords[2,]}))
    test_cont <- make_poly(pts_on_l, "test_cont")
    diff_reg1 <- gDifference(test_cont, conts[[j]])
    diff_reg2 <- gDifference(conts[[j]], test_cont)
    area_out[i, j] <- gArea(diff_reg1) + gArea(diff_reg2)
    print(c(i, j))
  }
}

#save(area_out, file = "/users/hdirector/desktop/area_out.rda")



#plot 
pdf("/users/hdirector/desktop/test_pts.pdf")
plot(coords_scale[[1]], type = "l", xlim = c(0, 1), ylim = c(0, 1), xlab = "", 
     ylab = "")
for (i in 1:n_obs) {
  points(coords_scale[[i]], type = "l", xlim = c(0, 1), ylim = c(0, 1))
}
plot(grid_pts, col = 'blue', pch = 20, cex = .5, add = T)
dev.off()
# rbPal <- colorRampPalette(c('red','blue'))
# cols <- rbPal(10)[as.numeric(cut(max_area_out,breaks = 10))]
# points(grid_pts@coords, col = cols, pch = 20)


#Group data into decades
#1:10 = 1980:1989, 11:20 = 1990:1999, 21:30 = 2000:2009, 31:38 = 2010:2017 
dec_names <- c("1980-1989", "1990-1999", "2000-2009", "2010-2017")
dec_ind <- matrix(nrow = 4, ncol = 2, data = c(1, 10, 11, 20, 21, 30, 31, 38),
                     byrow = TRUE)
n_case <- nrow(dec_ind)
pdf("/users/hdirector/desktop/approxPerform.pdf")
par(mfrow = c(2, 3), oma = c(0, 1, 0, 0), mar = c(0, 4, 1, 0))
for (i in 1:n_case) {
  inds_i <- dec_ind[i, 1]:dec_ind[i, 2]
  max_area_out <- apply(area_out[,inds_i], 1, max) 
  opt_ind <- which.min(max_area_out)
  opt_areas_i <- area_out[opt_ind, inds_i]
  C_hat <- matrix(data = grid_pts@coords[opt_ind,], ncol = 2)
  l_pts <- cbind(C_hat[1] + r*cos(theta), C_hat[2] + r*sin(theta))
  l <- apply(l_pts, 1, function(x){make_line(C_hat, x, "l")})
  
  #plot 3 cases (min, max, and, one in the middle)
  n_inds <- length(inds_i)
  for (j in c(1, floor(n_inds/2), n_inds)) {
    curr <- which(order(opt_areas_i) == j)
    cont_ind <- dec_ind[i, 1] + curr - 1
    on_l <- lapply(l, function(x){gIntersection(conts[[cont_ind]], x)})
    pts_on_l <- t(sapply(on_l, function(x){x@lines[[1]]@Lines[[1]]@coords[2,]}))
    test_cont <- make_poly(pts_on_l, "test_cont")
    if (j == 1) {
      plot(conts[[cont_ind]], main = sprintf("%i-th Ranked Area", j), 
           xlim = c(.1, .9), ylim = c(.1, .9), lwd = 1, ylab = dec_names[i])
    } else {
      plot(conts[[cont_ind]], main = sprintf("%i-th Ranked Area", j), 
           xlim = c(.1, .9), ylim = c(.1, .9), lwd = 1)
    }
    plot(test_cont, add = T, border = 'blue', lwd = 1)
    points(pts_on_l, col = 'blue', pch = 20)
    points(C_hat, col = 'red',pch = 20)
  }
}
dev.off()
