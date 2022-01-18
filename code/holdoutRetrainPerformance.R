# Plot the model performance in hold-one-out cross validation when
# when assuming the alignment based on gibbs sampling from all proteins

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(broom)
library(lvplot)

rm(list = ls())

ORACLE <- 'False'
CHAIN <- '100'
ITER <- '15'
RAND_SEED <- '382738375'
MWID <- '6'
EDGES <- 'cutAApos_1.0_0.05_edgeCut_1.0_0.05'
AMINO <- c('A','C','D','E','F','G','H','I','K','L',
           'M','N','P','Q','R','S','T','V','W','Y')
ADD_RESIDUES <- FALSE
CORE_POS <- paste0('A.',c('2','3','4','5','47','50','51','54','55'))
#if (ADD_RESIDUES) {
#  CORE_POS <- paste0('A.',c('2','3','4','5','6','47','50','51','54','55'))
#}

inDir <- paste0('../results/cisbp-chuAll/structFixed1_grpHoldout_multinomial_ORACLE',ORACLE,'Chain',CHAIN,
                'Iter',ITER,'scaled50/')
infile <- paste0(inDir,'pccTable_underS_holdOneOut.txt')
infile.noHO <- paste0(inDir,'pccTable_underS.txt')
if (ADD_RESIDUES) {
  outdir <- paste0(inDir,'modelFitPlots_holdOneOut_addResCS2012/')
} else {
  outdir <- paste0(inDir,'modelFitPlots_holdOneOut/')
}
dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

fitInfo <- fread(infile)
fitInfo$pos <- factor(fitInfo$pos)
fitInfo$pccAgree <- ifelse(fitInfo$pcc >= 0.5, TRUE, FALSE)
fitInfo$testType <- 'hold-one-out'
fitInfo$dset <- 
  ifelse(grepl('^[0-9]{1,3}_[ACGT]{3}_', fitInfo$prot), 'Chu2012', 'CIS-BP')

fitInfo.addRes <- fread(paste0(inDir,'pccTable_underS_holdOneOut_addResCS2012.txt'))
fitInfo.addRes$pos <- factor(fitInfo.addRes$pos)
fitInfo.addRes$pccAgree <- ifelse(fitInfo.addRes$pcc >= 0.5, TRUE, FALSE)
fitInfo.addRes$testType <- 'hoo-addRes-CS2012'
fitInfo.addRes$dset <- 
  ifelse(grepl('^[0-9]{1,3}_[ACGT]{3}_', fitInfo.addRes$prot), 'Chu2012', 'CIS-BP')

fitInfo.noHO <- fread(infile.noHO)
fitInfo.noHO$pos <- factor(fitInfo.noHO$pos)
fitInfo.noHO$pccAgree <- ifelse(fitInfo.noHO$pcc >= 0.5, TRUE, FALSE)
fitInfo.noHO$testType <- 'fit'
fitInfo.noHO$dset <- 
  ifelse(grepl('^[0-9]{1,3}_[ACGT]{3}_', fitInfo.noHO$prot), 'Chu2012', 'CIS-BP')

# Plot all the basic info about *this fit*
g <- ggplot(fitInfo, aes(x = pos, y = pcc)) + 
  geom_boxplot(outlier.colour = 'white') + 
  geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'pcc_boxplots.pdf'),height = 4, width = 5)

# Make cumulative plot with percentage of columns with PCC >= x
pcc <- c()
fracPCCge <- c()
colType <- c()
for (x in seq(-1,1, by = 0.01)) {
  pcc <- c(pcc, x)
  fracPCCge <- c(fracPCCge,nrow(fitInfo[pcc >= x])/nrow(fitInfo))
  colType <-  c(colType,'hold-one-out')
}
for (x in seq(-1,1, by = 0.01)) {
  pcc <- c(pcc, x)
  fracPCCge <- c(fracPCCge,nrow(fitInfo.noHO[pcc >= x])/nrow(fitInfo.noHO))
  colType <-  c(colType,'fit')
}
fracPCCge.tab <- data.table(pcc = pcc, fracPCCge = fracPCCge, colType = colType)
g <- ggplot(fracPCCge.tab, aes(x = pcc, y = fracPCCge)) + 
  geom_line(size = 1, aes(col = colType)) + 
  #geom_point(data = fracPCCge.tab[pcc == 0.5], shape = 2) +
  geom_vline(xintercept = 0.5, lty = 'dashed') +
  #geom_hline(yintercept = 0.95, lty = 'dashed') +
  scale_color_brewer("", palette = 'Dark2') + 
  labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
  theme_bw()
ggsave(plot = g, file = paste0(outdir, 'pcc_fracGE.pdf'),height = 4, width = 5)

g <- ggplot(fracPCCge.tab[colType == 'hold-one-out'], aes(x = pcc, y = fracPCCge)) + 
  geom_line(size = 1) + 
  #geom_point(data = fracPCCge.tab[pcc == 0.5], shape = 2) +
  geom_vline(xintercept = 0.5, lty = 'dashed') +
  #geom_hline(yintercept = 0.95, lty = 'dashed') +
  #scale_color_brewer("", palette = 'Dark2') + 
  labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
  theme_bw()
ggsave(plot = g, file = paste0(outdir, 'pcc_fracGE_noFit.pdf'),height = 4, width = 4)

g <- ggplot(fitInfo, aes(x = pos, y = rmse)) + 
  geom_boxplot(outlier.colour = 'white') + 
  geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4) + 
  #geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "SSE (predicted vs. actual)") +
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
g <- ggplot(fitInfo, aes(x = pos, y = pcc)) + 
  geom_boxplot(outlier.size = 0.1, aes(fill = dset)) + 
  #geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4, aes(color = dset)) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
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


# Compare *this fit* to the fit when not doing holdout
if (ADD_RESIDUES){
  fitInfo.both <- rbind(fitInfo, fitInfo.addRes) 
} else {
  fitInfo.both <- rbind(fitInfo, fitInfo.noHO)
}
g <- ggplot(fitInfo.both, aes(x = pos, y = pcc, fill = testType, grp = testType)) + 
  geom_boxplot(outlier.colour = 'white') + 
  scale_fill_brewer("",palette = 'Dark2') + 
  #scale_color_brewer(palette = 'Dark2') + 
  geom_point(size = 0.5, alpha = 0.4, shape = 21,
             position = position_jitterdodge(jitter.height = 0, jitter.width = 0.1)) + 
  #geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'cmpToFit_pcc_boxplots.pdf'),height = 4, width = 5)

# Boxen plots version
g <- ggplot(fitInfo.both, aes(x = pos, y = pcc, fill = testType, grp = testType)) + 
  geom_lv(color = 'black', outlier.size = 1) + 
  scale_fill_brewer("",palette = 'Dark2') + 
  #scale_color_brewer(palette = 'Dark2') + 
  #geom_point(size = 0.5, alpha = 0.4, shape = 21,
  #           position = position_jitterdodge(jitter.height = 0, jitter.width = 0.1)) + 
  #geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'cmpToFit_pcc_boxen_plots.pdf'),height = 4, width = 5)

g <- ggplot(fitInfo.both[testType == 'hold-one-out'], aes(x = pos, y = pcc)) +#, fill = testType, grp = testType)) + 
  geom_boxplot(outlier.colour = 'white') + 
  geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4, col = 'gray20') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'cmpToFit_pcc_boxplots_noFit.pdf'),height = 4, width = 4)

# Boxen plot version
g <- ggplot(fitInfo.both[testType == 'hold-one-out'], aes(x = pos, y = pcc)) +#, fill = testType, grp = testType)) + 
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  #scale_fill_brewer("", palette = "Dark2") +
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'cmpToFit_pcc_boxen_plots_noFit.pdf'),height = 4, width = 4)

g <- ggplot(fitInfo.both[testType == 'hold-one-out'], aes(x = pos, y = pcc)) +#, fill = testType, grp = testType)) + 
  geom_violin(color = 'black')+#, outlier.size = 1) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'cmpToFit_pcc_violin_plots_noFit.pdf'),height = 4, width = 4)

tmp <- fitInfo.both[,.(fracAgree = length(which(pccAgree == TRUE))/.N),
                    by = c('pos','testType')]
tmp$testType <- factor(tmp$testType)
g <- ggplot(tmp, aes(x = pos, y = fracAgree, fill = testType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer("",palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'cmpToFit_pcc_fracAgree.pdf'),height = 4, width = 5)

# Compare *this fit* to the fit when not doing holdout (just on Chu data)
g <- ggplot(fitInfo.both[dset == 'Chu2012'], aes(x = pos, y = pcc, fill = testType, grp = testType)) + 
  geom_boxplot(outlier.colour = 'white') + 
  scale_fill_brewer("",palette = 'Dark2') + 
  #scale_color_brewer(palette = 'Dark2') + 
  geom_point(size = 0.5, alpha = 0.4, shape = 21,
             position = position_jitterdodge(jitter.height = 0, jitter.width = 0.1)) + 
  #geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'cmpToFit_pcc_boxplots_Chu2012Only.pdf'),height = 4, width = 5)
tmp <- fitInfo.both[dset == 'Chu2012',.(fracAgree = length(which(pccAgree == TRUE))/.N),
                    by = c('pos','testType')]
tmp$testType <- factor(tmp$testType)
g <- ggplot(tmp, aes(x = pos, y = fracAgree, fill = testType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer("",palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'cmpToFit_pcc_fracAgree_Chu2012Only.pdf'),height = 4, width = 5)

# Compare *this fit* to the fit when not doing holdout (exclude Chu data)
g <- ggplot(fitInfo.both[dset != 'Chu2012'], aes(x = pos, y = pcc, fill = testType, grp = testType)) + 
  geom_boxplot(outlier.colour = 'white') + 
  scale_fill_brewer("",palette = 'Dark2') + 
  #scale_color_brewer(palette = 'Dark2') + 
  geom_point(size = 0.5, alpha = 0.4, shape = 21,
             position = position_jitterdodge(jitter.height = 0, jitter.width = 0.1)) + 
  #geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'cmpToFit_pcc_boxplots_noChu2012.pdf'),height = 4, width = 5)
tmp <- fitInfo.both[dset != 'Chu2012',.(fracAgree = length(which(pccAgree == TRUE))/.N),
                    by = c('pos','testType')]
tmp$testType <- factor(tmp$testType)
g <- ggplot(tmp, aes(x = pos, y = fracAgree, fill = testType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer("",palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'cmpToFit_pcc_fracAgree_noChu2012.pdf'),height = 4, width = 5)

# Are certain residues in core sequence positions good indicators of 
# correct vs. incorrect predictions?
aaTabs <- list()
aaTabs <- fitInfo[pos == 1,c('prot','coreSeq'),with=FALSE]
for (i in 1:nchar(aaTabs$coreSeq[1])) {
  aaTabs[[CORE_POS[i]]] <- factor(substr(aaTabs$coreSeq,i,i),levels = AMINO)
}
# Compute number of times each aa is observed in each position
aaTabs.melt <- melt(aaTabs, id.vars = 'prot', measure.vars = CORE_POS,
                    variable.name = 'position', value.name = 'aa')
g <- ggplot(aaTabs.melt, aes(x = aa)) + 
  geom_bar(color = 'black') + 
  facet_wrap(~position, scales = 'free_y') + 
  labs(x = 'Amino acid') +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'aa_perPos_allData.pdf'),
       height = 6, width = 8)
aaTabs.counts <- aaTabs.melt[,.(nObs = .N), by = c('position','aa')]
aaTabs.counts$ge5 <- factor(ifelse(aaTabs.counts$nObs >= 5, 1, 0), levels = c(0,1))
aaTabs.counts$ge10 <- factor(ifelse(aaTabs.counts$nObs >= 10, 1, 0), levels = c(0,1))

# For hold-one-out examples, are we more likely to get a column
# wrong if at least one of the core residues has been observed fewer 
# than 5 (or 10) times
X.nObs.ge5 <- c()    # One value per protein
X.nObs.ge10 <- c()   # One value per protein
setkey(aaTabs.melt, prot)
setkey(fitInfo, prot)
for (p in unique(aaTabs.melt$prot)) {
  tmp <- merge(aaTabs.melt[prot == p], aaTabs.counts, by = c('position','aa'))
  X.nObs.ge5 <- c(X.nObs.ge5, c(ifelse(all(tmp$ge5 == 1),1, 0)))
  X.nObs.ge10 <- c(X.nObs.ge10, c(ifelse(all(tmp$ge10 == 1),1, 0)))
}
lmInfo <- NULL
for (bp in 1:6) {
  y <- ifelse(fitInfo[pos == bp]$pccAgree,1,0)
  X.expIC.ge0.5 <- ifelse(fitInfo[pos == bp]$ic.exp > 0.5, 1, 0)
  tmp1 <- tidy(glm(y ~ X.nObs.ge5, family = 'binomial'))
  tmp2 <- tidy(glm(y ~ X.nObs.ge10, family = 'binomial'))
  tmp3 <- tidy(glm(y ~ X.expIC.ge0.5, family = 'binomial'))
  tmp <- rbind(rbind(tmp1, tmp2), tmp3)
  tmp$position <- bp
  lmInfo <- rbind(lmInfo, tmp)
}
rm(tmp, tmp1, tmp2)
lmInfo <- data.table(lmInfo)[term != '(Intercept)']
lmInfo$signif <- ifelse(lmInfo$p.value < 0.01, TRUE, FALSE)
lmInfo$signif.bh <- ifelse(p.adjust(lmInfo$p.value, 'BH') < 0.01, TRUE, FALSE)
lmInfo$position <- factor(lmInfo$position)
g <- ggplot(lmInfo, aes(x = position, y = estimate, color = term)) +
  geom_point(aes(shape = signif.bh)) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Binding site position", y = "Change in log-odds of column agreement") + 
  theme_classic()
ggsave(plot = g, 
       file = paste0(outdir, 'agreePCCvsNumObs_logRegression_fits.pdf'),
       height = 4, width = 7)

# Plot experimental IC distributions for correct vs incorrect columns
g <- ggplot(fitInfo, aes(x = pos, y = ic.exp, fill = pccAgree))+
  geom_boxplot(outlier.colour = 'white') + 
  scale_fill_brewer("Correct\nprediction",palette = 'Dark2') + 
  geom_point(size = 0.5, alpha = 0.4, shape = 21,
             position = position_jitterdodge(jitter.height = 0, jitter.width = 0.1)) + 
  labs(x = "Binding site position", y = "IC of experimental column") +
  theme_classic()
ggsave(plot = g, file = paste0(outdir, 'agreeVsIC_boxplots.pdf'),
       height = 4, width = 7)
