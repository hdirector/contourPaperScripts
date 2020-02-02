#load task_table
load("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/nObsnGenTaskTable.rda")
n_tasks <- nrow(task_table)
n_evals <- 50

#make table to store results
sum_table <- data.frame("n_obs" = rep(task_table$n_obs, each = 3),
                        "n_gen" = rep(task_table$n_gen, each = 3),
                        "shape_name"= rep(task_table$shape_name, each = 3),
                        "nom_cover" = rep(c(.8, .9, .95), n_tasks),
                        "med_cover" = rep(NA, 3*n_tasks),
                        "sd_cover "= rep(NA, 3*n_tasks),
                        "min_cover "= rep(NA, 3*n_tasks),
                        "max_cover "= rep(NA, 3*n_tasks))


for (i in c(1:2, 4:5, 7:8,10:11, 13:14, 16:17, 19:27)){#1:n_tasks) {
  n_obs_i <- task_table$n_obs[i]
  n_gen_i <- task_table$n_gen[i]
  shape_name_i <- task_table$shape_name[i]
  load(sprintf("/Users/hdirector/Dropbox/Contours/ContourPaperScripts/Simulations/sim_results/nObsnGen/nObsnGen_id%i_%s_obs%i_gen%i.rda",
               i, shape_name_i, n_obs_i, n_gen_i))
  cover <- res_cover$cover/n_evals
  curr_80 <- which(sum_table$n_obs == n_obs_i & sum_table$n_gen == n_gen_i &
                   sum_table$shape_name == shape_name_i & 
                   sum_table$nom_cover == .8)
  sum_table[curr_80, 5:8] <- c(median(cover[,1]), sd(cover[,1]), min(cover[,1]), 
                               max(cover[,1]))
  curr_90 <- which(sum_table$n_obs == n_obs_i & sum_table$n_gen == n_gen_i &
                   sum_table$shape_name == shape_name_i &
                   sum_table$nom_cover == .9)
  sum_table[curr_90, 5:8] <- c(median(cover[,2]), sd(cover[,2]), min(cover[,2]), 
                               max(cover[,2]))
  curr_95 <- which(sum_table$n_obs == n_obs_i & sum_table$n_gen == n_gen_i &
                     sum_table$shape_name == shape_name_i &
                     sum_table$nom_cover == .95)
  sum_table[curr_95, 5:8] <- c(median(cover[,3]), sd(cover[,3]), min(cover[,3]), 
                               max(cover[,3]))
}


#Number of observations
#sum_table$n_obs <- sum_table$n_obs
sum_table$nom_cover <- as.factor(sum_table$nom_cover)
library("tidyverse")
pdf("/Users/hdirector/Dropbox/Contours/Presentations/WGWinter2019/Figures/cover_samplesize.pdf",
    height = 4, width = 8)
ggplot(filter(sum_table, n_gen == 100),
       aes(x = jitter(n_obs, factor = .2), y = med_cover, group = shape_name, 
           col = nom_cover, shape = shape_name)) +
  geom_point(size = 3) +
  #ylim(.7, 1.05) + 
  ggtitle("Median Coverage vs. Sample Size") +
  geom_hline(yintercept = .8) +
  geom_hline(yintercept = .9) +
  geom_hline(yintercept = .95) +
  xlab("Number of observations (jittered)") +
  ylab("Coverage")
dev.off()
  
pdf("/Users/hdirector/Dropbox/Contours/Presentations/WGWinter2019/Figures/cover_ngen.pdf",
    height = 4, width = 8)
ggplot(filter(sum_table, n_obs == 20),
       aes(x = jitter(n_gen), y = med_cover, group = shape_name,
           col = nom_cover, shape = shape_name)) +
  geom_point(size = 3) +
  ylim(.7, 1.05) + 
  ggtitle("Median Coverage vs. Generated Contours") +
  geom_hline(yintercept = .8) +
  geom_hline(yintercept = .9) +
  geom_hline(yintercept = .95) +
  xlab("Number of Generated Contours (jittered)") +
  ylab("Coverage")
dev.off()

#Default table
temp <- filter(sum_table, (shape_name == "stop_sign"|shape_name == "tie"), n_obs == 20, 
       n_gen == 100)
library("xtable")
xtable(temp[,3:6], digits = 3)

sum_table[sum_table$n_obs == 20 & sum_table$n_gen == 100,]

  
  