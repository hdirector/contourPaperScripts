#Note: sims on cluster split into groups of 10 denoted
#a, b, c, d with seeds 101, 102, 103, and 104 respectively.
#Splitting done for speed up and memory reduction

rm(list = ls())
task_path <- "/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/misspecTaskTable.rda"
load(task_path)

n_tasks <- nrow(task_table)
n_evals <- 40

res_misspec <- data.frame("loc" = c(rep("rand", n_tasks/2), rep("fixed", n_tasks/2)),
                          "nom" = rep(c("0.8", "0.9", "0.95"), 2),
                          "mean_0_1" = rep(NA, n_tasks), "sd_0_1" = rep(NA, n_tasks),
                          "mean_2_3" = rep(NA, n_tasks), "sd_2_3" = rep(NA, n_tasks),
                          "mean_4_5" = rep(NA, n_tasks),"sd_4_5" = rep(NA, n_tasks))



for (i in 1:n_tasks) {
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/misspec/misspec_id%i_%s_nCurlMin%i_%s.rda",
               i, task_table[i,]$shape, task_table[i,]$n_curl_min, "a"))
  cover <- res$res_cover$cover
  for (j in c("b", "c", "d")) {
    load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/misspec/misspec_id%i_%s_nCurlMin%i_%s.rda",
          i, task_table[i,]$shape, task_table[i,]$n_curl_min, j))
    cover <- cover + res$res_cover$cover
  }
  row_ind <- ifelse(task_table[i,]$rand_loc == TRUE, 1, 4) 
  row_ind <- row_ind:(row_ind + 2)
  if (task_table[i,]$n_curl_min == 0) {
    res_misspec[row_ind,"mean_0_1"] <- apply(cover, 2, mean)/n_evals
    res_misspec[row_ind,"sd_0_1"] <- apply(cover, 2, sd)/n_evals
  } else if (task_table[i,]$n_curl_min == 2) {
    res_misspec[row_ind,"mean_2_3"] <- apply(cover, 2, mean)/n_evals
    res_misspec[row_ind,"sd_2_3"] <- apply(cover, 2, sd)/n_evals
  } else if (task_table[i,]$n_curl_min == 4) {
    res_misspec[row_ind,"mean_4_5"] <- apply(cover, 2, mean)/n_evals
    res_misspec[row_ind,"sd_4_5"] <- apply(cover, 2, sd)/n_evals
  }
  rm(cover)
}

xtable::xtable(res_misspec, digits = 2)


#coverage along transects
c1 <- 4; c2 <- 6
load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/misspec/misspec_id%i_%s_nCurlMin%i_%s.rda",
             c1, task_table[c1,]$shape, task_table[c1,]$n_curl_min, "a"))
cover_c1 <- res$res_cover$cover
for (j in c("b", "c", "d")) {
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/misspec/misspec_id%i_%s_nCurlMin%i_%s.rda",
               c1, task_table[c1,]$shape, task_table[c1,]$n_curl_min, j))
  cover_c1 <- cover_c1 + res$res_cover$cover
}

load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/misspec/misspec_id%i_%s_nCurlMin%i_%s.rda",
             c2, task_table[c2,]$shape, task_table[c2,]$n_curl_min, "a"))
cover_c2 <- res$res_cover$cover
for (j in c("b", "c", "d")) {
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/misspec/misspec_id%i_%s_nCurlMin%i_%s.rda",
               c2, task_table[c2,]$shape, task_table[c2,]$n_curl_min, j))
  cover_c2 <- cover_c2 + res$res_cover$cover
}
pdf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Figures/misspec_transect.pdf",
    height = 3, width = 5)
par(oma = rep(0, 4), mar = c(4, 4, 1, 1))
plot(cover_c1[,2]/n_evals, type= 'l', ylim = c(0, 1), xlab = "", ylab = '')
mtext(text = "Index of Test Line", side = 1, line = 2.5)
mtext(text = "Proportion Covered", side = 2, line = 2.5)
points(cover_c2[,2]/n_evals, type = 'l', col = 'blue')
abline(h = .9, lty = 2, col = 'red')
legend("bottom", ncol = 3, lty = c(1, 1, 2), col = c("black", "blue", "red"),
       legend = c("Random", "Fixed", "Expected"), cex = .8)
dev.off()

