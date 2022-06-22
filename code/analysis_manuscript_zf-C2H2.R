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
# Uses output from 'holdoutRetrainGLM_zfC2H2.py', where the GLM was retrained with
# each distinct set of binding amino acid residue combinations' protiens held out separately,
# to understand the expected de novo predictive performance of rCLAMPS in 
# a strict hold-out cross-validation setup.
rm(list = ls())

MWID <- '4'
AMINO <- c('A','C','D','E','F','G','H','I','K','L',
           'M','N','P','Q','R','S','T','V','W','Y')

inDir <- paste0('../my_results/zf-C2H2_250_50_seedFFSdiverse6/')
infile <- paste0(inDir,'pccTable_underS_holdOneOut.txt')
aliFile <- paste0(inDir,'registrationInfo.txt')
outdir <- paste0(inDir, 'plots/')
dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

makePCCgeTable <- function(fitInfo, icThresh) {
  # Trim the low IC edges from aligned experimental PWMs
  pwmTrim <- data.table(prot = unique(fitInfo$prot))
  pwmTrim$start <- sapply(pwmTrim$prot, function(x) {min(which(fitInfo[prot == x]$ic.exp >= icThresh))})
  pwmTrim$end <- sapply(pwmTrim$prot, function(x) {max(which(fitInfo[prot == x]$ic.exp >= icThresh))})
  fitInfo <- merge(fitInfo, pwmTrim, by = 'prot')
  fitInfo <- fitInfo[pos >= start & pos <= end]
  
  # Make cumulative table for percentage of columns with PCC >= x
  pcc <- c()
  fracPCCge <- c()
  fracPCCge.ffs <- c()
  fracPCCge.cbp <- c()
  colType <- c()
  for (x in seq(-1,1, by = 0.01)) {
    pcc <- c(pcc, x)
    fracPCCge <- c(fracPCCge,nrow(fitInfo[pcc >= x])/nrow(fitInfo))
    fracPCCge.ffs <- c(fracPCCge.ffs,nrow(fitInfo[pcc >= x & dset == 'ffs'])/nrow(fitInfo[dset == 'ffs']))
    fracPCCge.cbp <- c(fracPCCge.cbp,nrow(fitInfo[pcc >= x & dset == 'cisBP'])/nrow(fitInfo[dset == 'cisBP']))
    colType <-  c(colType,'hold-one-out')
  }
  fracPCCge.tab <- data.table(pcc = pcc, fracPCCge = fracPCCge, colType = colType, dset = 'all')
  tmp <- data.table(pcc = pcc, fracPCCge = fracPCCge.ffs, colType = colType, dset = 'ffs')
  fracPCCge.tab <- rbind(fracPCCge.tab, tmp)
  tmp <- data.table(pcc = pcc, fracPCCge = fracPCCge.cbp, colType = colType, dset = 'cisBP')
  fracPCCge.tab <- rbind(fracPCCge.tab, tmp)
  list("fitInfo" = fitInfo, "pccGEtab" = fracPCCge.tab)
}

# Fit when using strict holdout validation setup
fitInfo <- fread(infile)
fitInfo$domPos <- (fitInfo$pos-1)%%3 + 1
fitInfo$domPos <- factor(fitInfo$domPos)
fitInfo$pccAgree <- ifelse(fitInfo$pcc >= 0.5, TRUE, FALSE)
fitInfo$testType <- 'hold-one-out'
fitInfo$dset <- ifelse(grepl("^T", fitInfo$prot), "cisBP", "ffs")

trimRes.0 <- makePCCgeTable(fitInfo, icThresh = 0)
trimRes.25 <- makePCCgeTable(fitInfo, icThresh = 0.25)
trimRes.4 <- makePCCgeTable(fitInfo, icThresh = 0.4)


fitInfo.trim.0 <- trimRes.0[["fitInfo"]]
fracPCCge.tab.trim.0 <- trimRes.0[["pccGEtab"]]
fitInfo.trim.25 <- trimRes.25[["fitInfo"]]
fracPCCge.tab.trim.25 <- trimRes.25[["pccGEtab"]]
fitInfo.trim.4 <- trimRes.4[["fitInfo"]]
fracPCCge.tab.trim.4 <- trimRes.4[["pccGEtab"]]

# Figure 2 (left)
g <- ggplot(fracPCCge.tab.trim.0, aes(x = pcc, y = fracPCCge, col = dset)) + 
  geom_line(size = 1) + 
  #geom_point(data = fracPCCge.tab[pcc == 0.5], shape = 2) +
  geom_vline(xintercept = 0.5, lty = 'dashed') +
  #geom_hline(yintercept = 0.95, lty = 'dashed') +
  scale_color_brewer("", palette = 'Dark2') + 
  labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
  theme_bw()
ggsave(plot = g, file = paste0(outdir, 'Figure2_left.pdf'),height = 4, width = 5)

g <- ggplot(fracPCCge.tab.trim.25, aes(x = pcc, y = fracPCCge, col = dset)) + 
  geom_line(size = 1) + 
  #geom_point(data = fracPCCge.tab[pcc == 0.5], shape = 2) +
  geom_vline(xintercept = 0.5, lty = 'dashed') +
  #geom_hline(yintercept = 0.95, lty = 'dashed') +
  scale_color_brewer("", palette = 'Dark2') + 
  labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
  theme_bw()
ggsave(plot = g, file = paste0(outdir, 'Figure2_left_trimIC_0.25.pdf'),height = 4, width = 5)

g <- ggplot(fracPCCge.tab.trim.4, aes(x = pcc, y = fracPCCge, col = dset)) + 
  geom_line(size = 1) + 
  #geom_point(data = fracPCCge.tab[pcc == 0.5], shape = 2) +
  geom_vline(xintercept = 0.5, lty = 'dashed') +
  #geom_hline(yintercept = 0.95, lty = 'dashed') +
  scale_color_brewer("", palette = 'Dark2') + 
  labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
  theme_bw()
ggsave(plot = g, file = paste0(outdir, 'Figure2_left_trimIC_0.4.pdf'),height = 4, width = 5)

# Figure 2 (right)
g <- ggplot(fitInfo.trim.0, aes(x = domPos, y = pcc)) +
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'Figure2_right.pdf'),height = 4, width = 4)

g <- ggplot(fitInfo.trim.25, aes(x = domPos, y = pcc)) +
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'Figure2_right_trimIC_0.25.pdf'),height = 4, width = 4)

g <- ggplot(fitInfo.trim.4, aes(x = domPos, y = pcc)) +
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'Figure2_right_trimIC_0.4.pdf'),height = 4, width = 4)

# Alignment info 
aliInfo <- fread(aliFile)  # The alignment inferred by the procedure
aliFFS <- fread('../flyFactorSurvey/enuameh/enuameh_startPosInfo.txt')  # Optimal alignment with FFS
aliComp <- merge(aliFFS, aliInfo, by = 'prot')[prot %in% fitInfo$prot]
print(paste("We infer accurate registrations for", nrow(aliComp[start.x == start.y & rev.x == rev.y]),"out of",
            nrow(aliComp),"C2H2-ZFs with known registrations."))


