# A script for basic exploration of coefficients of learned by rCLAMPS

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)
rm(list = ls())

DOMAIN_TYPE <- 'homeodomain'  #Set to 'homeodomain' or 'zf-C2H2'
PLOT_TYPE <- 'coef' # set to 'coef' or 'offset'   # NOTE:  'offset' adds the intercept of the regression to each of the coeficient for that regression
PLOT_DIR <- paste0('../coefficientExplore/',DOMAIN_TYPE,'/',PLOT_TYPE,'_plots/')
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

x <- fread(paste0('../coefficientExplore/',DOMAIN_TYPE,'/coefTable.txt'))
contacts.all <- x[complete.cases(x)]
intercepts <- x[is.na(aapos)]
intercepts$intercept <- intercepts$coef
intercepts <- intercepts[,c('bpos', 'base', 'intercept'),with=FALSE]
contacts.all <- merge(contacts.all, intercepts, by = c('bpos','base'))
contacts.all$offset <- contacts.all$coef + contacts.all$intercept
contacts.all$contactPair <- paste(paste0('b',contacts.all$bpos), 
                                  paste0('a',contacts.all$aapos), sep = '.')
contacts.all$contactPair <- factor(contacts.all$contactPair)

for (cp in levels(contacts.all$contactPair)) {
  contacts <- contacts.all[contactPair == cp]
  
  # Create a 4 X 19 grid of coefficients and plot of relative offsets
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(32)
  g <- ggplot(droplevels(contacts), aes_string(x = 'aa', y = 'base', fill = PLOT_TYPE)) + 
    scale_fill_gradient2(low = muted("blue"), high = muted("red")) +
    geom_tile(color = 'black')
  ggsave(g, file=paste0(PLOT_DIR,PLOT_TYPE,'_',cp,'.pdf'), width = 8, height = 2)
}
