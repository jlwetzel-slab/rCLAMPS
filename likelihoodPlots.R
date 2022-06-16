# Examine the disctibution of best likelihood scores across Markov chains

library(data.table)
library(ggplot2)
library(RColorBrewer)

rm(list = ls())

IN_DIR <- './results/cisbp-chuAll/'

scores <- fread(paste0(IN_DIR, 'chainScores.txt'))
scores <- scores[order(-score)]
scores$rank <- 1:nrow(scores)

g <- ggplot(scores, aes(x = rank, y = score)) + 
  geom_point() +
  geom_line()
  

