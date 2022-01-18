library(data.table)
library(ggplot2)
library(RColorBrewer)

rm(list = ls())

ORACLE <- 'False'
CHAIN <- '100'
ITER <- '15'
RAND_SEED <- '382738375'
MWID <- '6'
EDGES <- 'cutAApos_1.0_0.05_edgeCut_1.0_0.05'

inDir <- paste0('../results/cisbp-chuAll/structFixed1_grpHoldout_multinomial_ORACLE',ORACLE,'Chain',CHAIN,
                'Iter',ITER,'scaled50/')
infile <- paste0(inDir,'pccTable_underS.txt')
outdir <- paste0(inDir,'modelFitPlots/')
dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

fitInfo <- fread(infile)
fitInfo$pos <- factor(fitInfo$pos)
fitInfo$pccAgree <- ifelse(fitInfo$pcc >= 0.5, TRUE, FALSE)


g <- ggplot(fitInfo, aes(x = pos, y = pcc)) + 
  geom_boxplot(outlier.colour = 'white') + 
  geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC Agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'pcc_boxplots.pdf'),height = 4, width = 5)

g <- ggplot(fitInfo, aes(x = pos, y = rmse)) + 
  geom_boxplot(outlier.colour = 'white') + 
  geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4) + 
  #geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'rmse_boxplots.pdf'),height = 4, width = 5)

g <- ggplot(fitInfo[pccAgree == TRUE], aes(x = pos)) + 
  geom_bar(color = 'black', width = 0.6) + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'pcc_numAgree.pdf'),height = 4, width = 5)

tmp <- fitInfo[,.(fracAgree = length(which(pccAgree == TRUE))/.N),by = c('pos')]
g <- ggplot(tmp, aes(x = pos, y = fracAgree)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity') + 
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer(palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'pcc_fracAgree.pdf'),height = 4, width = 5)

# How do things look for the two different datasets?
fitInfo$dset <- 
  ifelse(grepl('^[0-9]{1,3}_[ACGT]{3}_', fitInfo$prot), 'Chu2012', 'CIS-BP')
g <- ggplot(fitInfo, aes(x = pos, y = pcc)) + 
  geom_boxplot(outlier.size = 0.1, aes(fill = dset)) + 
  #geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4, aes(color = dset)) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC Agreement (predicted vs. actual)") +
  scale_fill_brewer(palette = 'Dark2') + 
  #scale_color_brewer(palette = 'Dark2') +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'pcc_boxplots_dsetColors.pdf'),height = 4, width = 5)

tmp <- fitInfo[,.(fracAgree = length(which(pccAgree == TRUE))/.N),
               by = c('pos','dset')]
g <- ggplot(tmp, aes(x = pos, y = fracAgree, fill = dset)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer(palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'pcc_fracAgree_dsetColors.pdf'),height = 4, width = 5)
