# Compare results from optimal output of Gibbs sampler compared 
# against STAMP

library(data.table)
library(ggplot2)
library(RColorBrewer)

rm(list = ls())

ORACLE <- FALSE
CHAIN <- '100'
ITER <- '15'
DSET_LAB <- c('STAMP', 'rCLAMP')

plotDir <- paste0('../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLE',
                  ORACLE,'Chain',CHAIN,'Iter',ITER,'/flyOrthologAlignmentSummary/')
dir.create(plotDir,showWarnings = FALSE,recursive = TRUE)
niceCols <- RColorBrewer::brewer.pal(8, "Dark2")#[1:3]

# Plot number validated as correct per visual inspection of alignment
sgrps <- c('Abd-B','Antp','Bar','Bcd','En','Iroquois','Ladybird',
           'NK-1','NK-2','Six','TGIF-Exd')
specGrps <- rep(sgrps,3)
valGrps <- c(rep(DSET_LAB[1],length(sgrps)),
             rep(DSET_LAB[2],length(sgrps)),
             rep('possible',length(sgrps)))
vals <- c(14, 44, 15, 5, 73, 0, 5, 14, 6, 5, 6,
          35, 43, 17, 5, 78, 0, 5, 20, 14, 5, 6,
          36, 45, 18, 5, 78, 3, 5, 20, 14, 5, 6)

mappingStats <- data.table(specGrp = specGrps, valType = valGrps, value = vals)
g <- ggplot(droplevels(mappingStats[valGrps != 'possible']), 
            aes(x = specGrp, y = value)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'black',
           width = 0.7, aes(fill = valType)) + 
  geom_point(shape = 1, data = mappingStats[valGrps == 'possible']) + 
  scale_fill_manual("Method", values = niceCols) + 
  labs(x = 'Specificity group', y = '# Correct mappings') + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'top', legend.title = element_blank())
ggsave(plot = g, file = paste0(plotDir, 'mappingsCorrect.pdf'),
       height = 3, width = 3.5)

g <- ggplot(droplevels(mappingStats[valGrps != 'possible']), 
            aes(x = specGrp, y = value)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'black',
           width = 0.7, aes(fill = valType)) + 
  geom_point(shape = 1, data = mappingStats[valGrps == 'possible']) + 
  scale_fill_manual("Method", values = niceCols) + 
  labs(x = 'Specificity group', y = '# Correct mappings') + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none', legend.title = element_blank())
ggsave(plot = g, file = paste0(plotDir, 'mappingsCorrect_noLegend.pdf'),
       height = 3, width = 3.5)

mappingStats.summ <- 
  data.table(fracCorrect = 
               c(sum(mappingStats[valType == DSET_LAB[1]]$value)/
                   sum(mappingStats[valType == 'possible']$value),
                 sum(mappingStats[valType == DSET_LAB[2]]$value)/
                   sum(mappingStats[valType == 'possible']$value)),
             dset = DSET_LAB)
g <- ggplot(mappingStats.summ, aes(x = dset, y = fracCorrect)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'black',
           width = 0.7, aes(fill = dset)) + 
  scale_fill_manual("Method", values = niceCols) + 
  scale_y_continuous(limits = c(0,1)) +
  labs(x = 'Method', y = 'Fraction of correct mappings') + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),legend.position = 'none')
ggsave(plot = g, file = paste0(plotDir, 'mappingsCorrect_summary.pdf'),
       height = 3, width = 2)