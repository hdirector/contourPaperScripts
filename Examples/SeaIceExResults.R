rm(list = ls())
years <- 2010:2017
p <- 66
cover_all <- matrix(nrow = p, ncol = 3, data = 0)
for (y in years) {
  load(sprintf("/Users/hdirector/Desktop/cover_year%i.rda", y))
  cover_all <- cover_all + cover
}

df <- data.frame("nominal" = c(.8, .9, .95),
                "mean" = apply(cover_all/8, 2, mean),
                "sd" = apply(cover_all/8, 2, sd))
library("xtable")
xtable(df, digits = 3)
