years <- 2008:2017
n_years <- length(years)
p <- 71
cover_all <- matrix(nrow = p, ncol = 3, data = 0)
for (year in years) {
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Examples/ex_results/cover_year%i.rda", year))
  cover_all <- cover_all + cover
}


df <- data.frame("nominal" = c(.8, .9, .95),
                "mean" = apply(cover_all/n_years, 2, mean),
                "sd" = apply(cover_all/n_years, 2, sd))
library("xtable")
xtable(df, digits = 3)

