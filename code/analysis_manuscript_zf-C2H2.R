# Reproduction of Figures in the zf-C2H2 figures for article

# Plot the model performance in hold-one-out cross validation when
# when assuming the alignment based on gibbs sampling from all proteins

library(plyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(broom)
library(lvplot)

###########################
# Figure 2
###########################
# Uses output from 'holdoutRetrainGLM.py', where the GLM was retrained with
# each distinct set of binding amino acid residue combinations' protiens held out separately,
# to understand the expected de novo predictive performance of rCLAMPS in 
# a strict hold-out cross-validation setup.
rm(list = ls())

MWID <- '4'
AMINO <- c('A','C','D','E','F','G','H','I','K','L',
           'M','N','P','Q','R','S','T','V','W','Y')

inDir <- paste0('../my_results/zf-C2H2_100_25_seedFFSall/')
infile <- paste0(inDir,'pccTable_underS_holdOneOut.txt')
outdir <- paste0(inDir, 'plots/')
dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

# Fit when using strict holdout validation setup
fitInfo <- fread(infile)
fitInfo$domPos <- (fitInfo$pos-1)%%3
fitInfo$pos <- factor(fitInfo$pos)
fitInfo$domPos <- factor(fitInfo$domPos)
fitInfo$pccAgree <- ifelse(fitInfo$pcc >= 0.5, TRUE, FALSE)
fitInfo$testType <- 'hold-one-out'

# Make cumulative plot with percentage of columns with PCC >= x
pcc <- c()
fracPCCge <- c()
colType <- c()
for (x in seq(-1,1, by = 0.01)) {
  pcc <- c(pcc, x)
  fracPCCge <- c(fracPCCge,nrow(fitInfo[pcc >= x])/nrow(fitInfo))
  colType <-  c(colType,'hold-one-out')
}
fracPCCge.tab <- data.table(pcc = pcc, fracPCCge = fracPCCge, colType = colType)

# Figure 2 (left)
g <- ggplot(fracPCCge.tab[colType == 'hold-one-out'], aes(x = pcc, y = fracPCCge)) + 
  geom_line(size = 1) + 
  #geom_point(data = fracPCCge.tab[pcc == 0.5], shape = 2) +
  geom_vline(xintercept = 0.5, lty = 'dashed') +
  #geom_hline(yintercept = 0.95, lty = 'dashed') +
  #scale_color_brewer("", palette = 'Dark2') + 
  labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
  theme_bw()
ggsave(plot = g, file = paste0(outdir, 'Figure2_left.pdf'),height = 4, width = 4)

# Figure 2 (right)
g <- ggplot(fitInfo[testType == 'hold-one-out'], aes(x = domPos, y = pcc)) +
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'Figure2_right.pdf'),height = 4, width = 4)

###############################
## 
###############################