# A script for basic exploration of coefficients of learned by rCLAMPS

library(data.table)
library(ggplot2)
library(RColorBrewer)
rm(list = ls())

DOMAIN_TYPE = 'zf-C2H2'

x <- fread(paste0('../examplePredictions/',DOMAIN_TYPE,'/coefTable.txt'))
#x <- x[complete.cases(x)]

# Create a 4 X 19 grid of coefficients and plot relative to 0
ref.base = 'T'; bp <- 1; ap <- 6
contacts <- x[bpos == bp & aapos == ap]
#contacts$coef <- sapply(1:nrow(contacts), function(i) {
#  ref.delta <- contacts[base == ref.base & aa == contacts$aa[i]]$coef
#  contacts$coef[i] - ref.delta
#})
#contacts <- contacts[base != ref.base]

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
g <- ggplot(droplevels(contacts), aes(x = aa, y = base, fill = coef)) + 
  scale_fill_gradient2() +
  geom_tile(color = 'black')
