# #Figure showing that data are approximately star-shaped
# pdf("/users/hdirector/desktop/approxPerform.pdf")
# par(mfrow = c(2, 3), oma = c(0, 1, 0, 0), mar = c(0, 4, 1, 0))
# for (i in 1:n_case) {
#    l <- make_l(C = C_hat[i,], theta)
#   #plot 3 cases (min, max, and, one in the middle)
#   n_inds <- length(inds_i)
#   for (j in c(1, floor(n_inds/2), n_inds)) {
#     curr <- which(order(opt_areas_i) == j)
#     cont_ind <- dec_ind[i, 1] + curr - 1
#     on_l <- lapply(l, function(x){gIntersection(conts[[cont_ind]], x)})
#     pts_on_l <- t(sapply(on_l, function(x){x@lines[[1]]@Lines[[1]]@coords[2,]}))
#     test_cont <- make_poly(pts_on_l, "test_cont")
#     if (j == 1) {
#       plot(conts[[cont_ind]], main = sprintf("%i-th Ranked Area", j), 
#            xlim = c(.1, .9), ylim = c(.1, .9), lwd = 1, ylab = dec_names[i])
#     } else {
#       plot(conts[[cont_ind]], main = sprintf("%i-th Ranked Area", j), 
#            xlim = c(.1, .9), ylim = c(.1, .9), lwd = 1)
#     }
#     plot(test_cont, add = T, border = 'blue', lwd = 1)
#     points(pts_on_l, col = 'blue', pch = 20)
#     points(C_hat[i, 1], C_hat[i, 2], col = 'red',pch = 20)
#   }
# }
# dev.off()