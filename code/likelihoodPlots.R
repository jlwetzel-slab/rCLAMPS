# Examine the disctibution of best likelihood scores across Markov chains

library(data.table)
library(ggplot2)
library(RColorBrewer)

rm(list = ls())

IN_DIR <- '../my_results/allHomeodomainProts/'
#IN_DIR <- '../my_results/zf-C2H2_250_50_seedFFSdiverse6/'

scores.rand <- fread(paste0(IN_DIR, 'randomModelScores.txt'))
#scores.rand$score <- log(scores.rand$score)/log(2)
scores <- fread(paste0(IN_DIR, 'chainScores.txt'))
scores <- scores[order(-score)]
scores$score <- scores$score - mean(scores.rand$score)
scores$rank <- 1:nrow(scores)

g <- ggplot(scores, aes(x = rank, y = score)) + 
  geom_point() +
  labs(y = "Improvement in log likelihood above random",
       x = "Rank") +
  theme_bw()
ggsave(g, file = paste0(IN_DIR, 'modelImprovementAboveRandom.pdf'),
       height = 4, width = 6)
  

