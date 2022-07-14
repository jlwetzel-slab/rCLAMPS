# Reproduction of Figures in the article

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

inDir <- paste0('../results/cisbp-chuAll/structFixed1_grpHoldout_multinomial_ORACLE',ORACLE,'Chain',CHAIN,
                'Iter',ITER,'scaled50/')
infile <- paste0(inDir,'pccTable_underS_holdOneOut.txt')
infile.noHO <- paste0(inDir,'pccTable_underS.txt')
outdir <- '../analysis_manuscript_plots/'
dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

# Fit when using strict holdout validation setup
fitInfo <- fread(infile)
fitInfo$pos <- factor(fitInfo$pos)
fitInfo$pccAgree <- ifelse(fitInfo$pcc >= 0.5, TRUE, FALSE)
fitInfo$testType <- 'hold-one-out'
fitInfo$dset <- 
  ifelse(grepl('^[0-9]{1,3}_[ACGT]{3}_', fitInfo$prot), 'Chu2012', 'CIS-BP')

# Fit without holdout
fitInfo.noHO <- fread(infile.noHO)
fitInfo.noHO$pos <- factor(fitInfo.noHO$pos)
fitInfo.noHO$pccAgree <- ifelse(fitInfo.noHO$pcc >= 0.5, TRUE, FALSE)
fitInfo.noHO$testType <- 'fit'
fitInfo.noHO$dset <- 
  ifelse(grepl('^[0-9]{1,3}_[ACGT]{3}_', fitInfo.noHO$prot), 'Chu2012', 'CIS-BP')

# Read in the info for the fit 
fitInfo <- fread(infile)
fitInfo$pos <- factor(fitInfo$pos)
fitInfo$pccAgree <- ifelse(fitInfo$pcc >= 0.5, TRUE, FALSE)
fitInfo$testType <- 'hold-one-out'
fitInfo$dset <- 
  ifelse(grepl('^[0-9]{1,3}_[ACGT]{3}_', fitInfo$prot), 'Chu2012', 'CIS-BP')

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

# Figure 2 (left)
g <- ggplot(fracPCCge.tab[colType == 'hold-one-out'], aes(x = pcc, y = fracPCCge)) + 
  geom_line(size = 1) + 
  #geom_point(data = fracPCCge.tab[pcc == 0.5], shape = 2) +
  geom_vline(xintercept = 0.5, lty = 'dashed') +
  #geom_hline(yintercept = 0.95, lty = 'dashed') +
  #scale_color_brewer("", palette = 'Dark2') + 
  labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
  theme_bw() + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
ggsave(plot = g, file = paste0(outdir, 'Figure2_left.pdf'),height = 4, width = 4)

# Figure 2 (right)
fitInfo.both <- rbind(fitInfo, fitInfo.noHO)
g <- ggplot(fitInfo.both[testType == 'hold-one-out'], aes(x = pos, y = pcc)) +
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
ggsave(plot = g, file = paste0(outdir, 'Figure2_right.pdf'),height = 4, width = 4)

###############################
## 
###############################




################
# Overall comparison to rf_extant and rf_joint, previously exising random forest methods 
# plus Figure 3 
################

# Fig 3 uses output from 'hd1_transfer_predictions.py', where hybrid approach was applied
# to wildtype (training) and neighboring (mutant/test) proteins understand value of 
# per-base-position transfer of specificity information vs. de novo prediction

rm(list = ls())
outdir <- '../analysis_manuscript_plots/'
RES_DIR <- '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/transfer_test/fracID_fullDBD/'
MATCHTABS_FILE <- '../hd1-explore/0_splitChu2012-trainTest/0_matchTabs.txt'
HD1_FILE <- '../hd1-explore/0_splitChu2012-trainTest/0_hd1cores_test_train.txt'
TRAIN_FIT_TAB <- paste0('../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/pccTable_underS.txt')
BARRERA_TAB <- '../barrera2016_SuppTable_S6_combined.csv'
RM_LOW_IC <- FALSE
RM_LOW_OBS_PREDS <- FALSE
CORE_POS <- paste0('A.',c('2','3','4','5','47','50','51','54','55'))
INCLUDE_RF_ALI <- FALSE#TRUE#

# Read in the hd1 info between test and train proteins
hd1Info <- fread(HD1_FILE)

# Read in the fit info from the training data
trainFit <- fread(TRAIN_FIT_TAB)

# Read in table S6 from Barrera 2016
s6tab <- data.table(read.csv(BARRERA_TAB))
s6tab$prot <- paste(s6tab$prot, s6tab$sub, sep = '_')

# Read in the prediction vs test pwm alignment results
resFiles <- c('pccTable_test_predOnly.txt', 'pccTable_test_transfer.txt')
resFileLabs <- c('model', 'model+nn')
resFiles_rf <- c('pccTable_test_pred.txt')
resFileLabs_rf <- c('rf_preAligned')
stormoFiles <- c('../stormo_predictor_outputs/extant/pccTable_test.txt', 
                 '../stormo_predictor_outputs/joint/pccTable_test.txt')
stormoLabs <- c('rf_extant','rf_joint')
res <- NULL
for (i in 1:length(resFiles)) {
  tmp <- fread(paste0(RES_DIR, resFiles[i]))
  tmp$predType <- resFileLabs[i]
  res <- rbind(res, tmp)
}
for (i in 1:length(stormoFiles)) {
  tmp <- fread(stormoFiles[i])
  tmp$predType <- stormoLabs[i]
  res <- rbind(res, tmp)
}
if (INCLUDE_RF_ALI) {
  for (i in 1:length(resFiles_rf)) {
    #i <- 1
    tmp <- fread(paste0(RES_DIR_RF, resFiles_rf[i]))
    tmp$predType <- resFileLabs_rf[i]
    res <- rbind(res, tmp)
  }
}
rm(tmp); res$pos <- factor(res$pos)
res$predType <- factor(res$predType, levels = c(resFileLabs, stormoLabs,resFileLabs_rf))
test.me <- unique(res[predType %in% c('model','model+transfer')]$prot)
test.stormo <- unique(res[predType %in% c('rf_extant','model+transfer')]$prot)
test.both <- intersect(test.me, test.stormo)

# Subset only to proteins where all could be predicted by both methods
res <- res[prot %in% test.both]

# Read in the match state and training/testing set info
matchTab <- fread(MATCHTABS_FILE)
matchTab.train <- matchTab[prot %in% unique(trainFit$prot) & grp == 'train']
matchTab.test <- matchTab[prot %in% unique(res$prot) & grp == 'test']

# Remove results for proteins where affinity was lost and few 
# k-mers were present to create PWMs
barrera.drop <- unique(s6tab[aff.change == '-' & n.8mers.alt < 10]$prot)
res <- res[!(prot %in% barrera.drop)]
res <- merge(res, matchTab[,c('prot','dset','coreSeq'),with=FALSE])

# Info for removing proteins used for training rf_extant or rf_joint
motifInfo.cisbp <- fread('~/research/cisBP/hboxOnlyInfo.txt')
trPubs.rfExtant <- c('Berger08','Zhu09','FlyFactorSurvey')
trSet.rfExtant <- unique(motifInfo.cisbp[MSource_Identifier %in% trPubs.rfExtant]$TF_Name)
#trSet.rfExtant.cores <- unique(res[predType != 'rf_joint' & !(prot %in% trSet.rfExtant))
trSet.rfJoint <- 
  union(unique(motifInfo.cisbp[MSource_Identifier %in% trPubs.rfExtant]$TF_Name),
        unique(matchTab[dset == 'Chu2012']$prot))

# Overall comparison with no proteins from test set from either method included
res.noOlap.rfExtant <- res[predType != 'rf_joint' & !(prot %in% trSet.rfExtant)]
res.noOlap.rfJoint <- res[predType != 'rf_extant' & !(prot %in% trSet.rfJoint)]
fracPCCge.list <- list()
dset.name <- c('rf_extant','rf_joint')
dset <- list(res.noOlap.rfExtant,res.noOlap.rfJoint)
for (i in 1:length(dset)) {
  tmp <- dset[[i]]
  pcc <- c()
  fracPCCge <- c()
  colType <- c()
  for (pt in unique(tmp$predType)) {
    for (x in seq(-1,1, by = 0.01)) {
      pcc <- c(pcc, x)
      fracPCCge <- c(fracPCCge,nrow(tmp[predType == pt & pcc >= x])/nrow(tmp[predType == pt]))
      colType <-  c(colType,pt)
    }
  }
  x <- data.table(pcc = pcc, fracPCCge = fracPCCge, colType = colType)
  fracPCCge.list[[dset.name[i]]] <- x
  # Exclude model+nn
  g <- ggplot(fracPCCge.list[[dset.name[i]]][colType != 'model+nn'], aes(x = pcc, y = fracPCCge)) + 
    geom_line(size = 0.75, aes(col = colType)) + 
    geom_vline(xintercept = 0.5, lty = 'dashed') +
    scale_color_brewer("", palette = 'Dark2') + 
    labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
    theme_bw()+ 
    theme(legend.position = 'top')
  #ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracGE_',dset.name[i],'.pdf'),height = 4, width = 6)
  ggsave(plot = g, file = paste0(outdir, 'pcc_fracGE_',dset.name[i],'_excludeTransfer.pdf'),height = 4.5, width = 4)
}

## Figure 3 - plot the overall comparison of transfer vs OUR model
res.fig3 <- merge(res,
                  res[predType == 'model+nn',.(transfer = nnbrs > 0),by=c('prot','pos')],
                  by = c('prot','pos'))
res.fig3$transfer <- ifelse(res.fig3$transfer, "transfer", "de novo")
nearMuts <- unique(res.fig3[transfer == 'transfer']$prot)
res.fig3 <- res.fig3[prot %in% nearMuts]
res.fig3$agree <- ifelse(res.fig3$pcc >= 0.5, TRUE, FALSE)
res.fig3$predType <- plyr::mapvalues(res.fig3$predType,
                                     from = c('model', 'model+nn','rf_extant','rf_joint'),
                                     to = c('model', 'hybrid','rf_extant','rf_joint'))
res.fig3$predType <- factor(res.fig3$predType,
                            levels = c('hybrid','model','rf_extant','rf_joint'))
res.fig3$predType <- plyr::revalue(res.fig3$predType, 
                                   replace = c('hybrid'='hybrid','model'='rCLAMPS',
                                               'rf_extant' = 'rf_extant','rf_joint'='rfJoint'))

# Figure 3 (left)
g <- ggplot(res.fig3[predType == 'hybrid'], aes(x = transfer, y = pcc)) + 
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  scale_y_continuous(limits = c(-1,1)) +
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Column type", y = "PCC between predicted and actual") +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
ggsave(plot = g, file = paste0(outdir, 'Figure3_left.pdf'),
       height = 4, width = 3)

# Figure 3 (right)
g <- ggplot(res.fig3[!(prot %in% trSet.rfExtant) &
                       predType %in% c('rCLAMPS','hybrid','rf_extant') & transfer == 'transfer'], 
            aes(x = predType, y = pcc)) + 
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  scale_y_continuous(limits = c(-1,1)) +
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Method", y = "PCC between predicted and actual") +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 11))
ggsave(plot = g, file = paste0(outdir, 'Figure3_right.pdf'),
       height = 4, width = 3)

################
# Figure 4
################
rm(list = ls())

ORACLE <- FALSE
CHAIN <- '100'
ITER <- '15'
DSET_LAB <- c('STAMP', 'rCLAMPS')

inDir <- paste0('../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLE',
                ORACLE,'Chain',CHAIN,'Iter',ITER,'/flyOrthologAlignmentSummary/')
dir.create(inDir,showWarnings = FALSE,recursive = TRUE)
niceCols <- RColorBrewer::brewer.pal(8, "Dark2")#[1:3]
outdir <- '../analysis_manuscript_plots/'

# Plot number validated as correct per visual inspection of alignment
sgrps <- c('Abd-B','Antp','Bar','Bcd','En','Iroquois','Ladybird',
           'NK-1','NK-2','Six','TGIF-Exd')
specGrps <- rep(sgrps,3)
valGrps <- c(rep(DSET_LAB[1],length(sgrps)),
             rep(DSET_LAB[2],length(sgrps)),
             rep('possible',length(sgrps)))

# Input from inspection of motifs with identical contacting amino acid 
# sequence identical to fly proteins (ref Noyes 2008b) and comparison to 
# experimental alignment (see ref. Noyes2008b, Figure2c and Supplemental Data).
# Motifs aligned by each method located in:
# '../STAMP/alignedLogos_NoyesGrps/logos_grpByNoyes08_hd0/' (STAMP)
# Aligned visual motifs for rCLAMPS were produced by gibbsAlign_PWM.py and located in:
# '../homeodomain/gibbsAlign_output_grpIDcore/cisbp-chu-scaled50/structFixed1_grpHoldout_glm_multinomial/cisbp_rSeed_382738375/mWid_6/nChains_100/cutAApos_1.0_0.05_edgeCut_1.0_0.05/logos_grpByNoyes08_hd0/' (rCLAMPS)
# STAMP was provided input located in '../STAMP/inputFile.txt' then its STAMP's output alignment was processed 
# using baselineAlignment.py to visualize the aligned motifs.
vals <- c(14, 44, 15, 5, 73, 0, 5, 14, 6, 5, 6,   # STAMP
          35, 43, 17, 5, 78, 0, 5, 20, 14, 5, 6,  # rCLAMPS
          36, 45, 18, 5, 78, 3, 5, 20, 14, 5, 6)  # Maximum

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
ggsave(plot = g, file = paste0(inDir, 'mappingsCorrect.pdf'),
       height = 3, width = 3.5)

g <- ggplot(droplevels(mappingStats[valGrps != 'possible']), 
            aes(x = specGrp, y = value)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'black',
           width = 0.7, aes(fill = valType)) + 
  geom_point(shape = 1, data = mappingStats[valGrps == 'possible']) + 
  scale_fill_manual("Method", values = niceCols) + 
  labs(x = 'Specificity group', y = '# Correct mappings') + 
  theme_classic() + 
  theme(axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        legend.title = element_blank(),legend.position = 'none')
ggsave(plot = g, file = paste0(outdir, 'Figure4_right.pdf'),
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
  theme(axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),legend.position = 'none')
ggsave(plot = g, file = paste0(outdir, 'Figure4_left.pdf'),
       height = 3, width = 2)

