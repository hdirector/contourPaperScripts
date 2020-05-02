#Assess coverage under particular conditions
rm(list = ls())
set.seed(103)
library("IceCast")
library("ContouR")
library("raster")
library("coda")
library("fields")
library("viridis")

#data about original grid
xmn = -3850; xmx = 3750; ymn = -5350; ymx = 5850
box_size <- 25
demo_year <- 2017

#load data
all <- IceCast::read_monthly_BS(demo_year, demo_year, 
                                "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/bootstrapV3_1/",
                                version = 3.1)
sep <- all[,9,,]

#find sea ice covered polygon and contour points
ice <- disaggregate(get_region(sep, dat_type = "bootstrap", level = 15))
ice_main <- ice[which.max(gArea(ice, byid = TRUE))]
ice_main <- disaggregate(rm_holes(keep_poly(gIntersection(reg_info$regions[[1]], ice_main))))
ice_main <- ice_main[which.max(gArea(ice_main, byid = T))]
ice_main_coords <- ice_main@polygons[[1]]@Polygons[[1]]@coords

#set region to consider
bb <-  matrix(ncol = 2, byrow = TRUE, data = c(-2300, -1000,
                                               -2300, 2100,
                                               1500, 2100, 
                                               1500, -1000))
bb <- SpatialPolygons(list(Polygons(list(Polygon(bb)), "bb")))


#plot motivating example figure
pdf(file = "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/figures/motivating_ex.pdf", 
    height = 3.5, width = 3.6)
par(oma = rep(0, 4), mar = c(2, 0, 0, 0))
plot(bb, col = 'blue')
plot(gDifference(bb, reg_info$regions[[1]]), col = 'beige', border = 'beige',
     add = T)
plot(gIntersection(land, bb), add = T, col = 'grey', border = "grey")
plot(gIntersection(ice, bb), add = T, col = 'white', border = 'white')

plot(bb, add = T)
plot(ice_main, border = 'red', add = TRUE, lwd = 1.5)
legend(-2300, -1100, cex = .65, ncol = 3,
       legend = c("contour", "ice", "no ice", "land", "outside region"),
       pch = c(NA, rep(22, 4)), lty = c(1, rep(NA, 4)),
       col = c("red", rep(NA, 4)),
       lwd = c(3, rep(NA, 4)), pt.cex = 1,
       text.font = 1, fill = c("white", "white", "blue", "grey", "beige"),
       xpd = NA,  bg = "white", border = c("white", rep("black",4)))
dev.off()

