library(data.table)
library(ggplot2)
library(RColorBrewer)
library(lvplot)
rm(list = ls())

RES_DIR <- '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/transfer_test/fracID_fullDBD/'
RES_DIR_RF <- '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/randomForest_regression_scaled50_optuna/testFit/'
#RES_DIR_RF <- '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/randomForest_regression_scaled50/structEdges/testFit/'
#RES_DIR <- '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15/transfer_test/mean/'
MATCHTABS_FILE <- '../hd1-explore/0_splitChu2012-trainTest/0_matchTabs.txt'
HD1_FILE <- '../hd1-explore/0_splitChu2012-trainTest/0_hd1cores_test_train.txt'
TRAIN_FIT_TAB <- paste0('../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/pccTable_underS.txt')
BARRERA_TAB <- '../barrera2016_SuppTable_S6_combined.csv'
RM_LOW_IC <- FALSE
RM_LOW_OBS_PREDS <- FALSE
CORE_POS <- paste0('A.',c('2','3','4','5','47','50','51','54','55'))
INCLUDE_RF_ALI <- FALSE#TRUE#

# Read in the table of model coefficients
mcoefs = fread(paste0(RES_DIR, 'modelCoefs.txt'))

# Read in the fit info from the training data
trainFit <- fread(TRAIN_FIT_TAB)

# Read in the hd1 info between test and train proteins
hd1Info <- fread(HD1_FILE)

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
tmp <- melt(matchTab.train, id.vars = c('prot','grp'), measure.vars = CORE_POS,
            variable.name = 'position', value.name = 'aa')
aaCounts.train <- tmp[,.(count = .N), by = c('position','aa')]
aaCounts.train$grp <- 'train'
tmp <- melt(matchTab.test, id.vars = c('prot','grp'), measure.vars = CORE_POS,
            variable.name = 'position', value.name = 'aa')
aaCounts.test <- tmp[,.(count = .N), by = c('position','aa')]
aaCounts.test$grp <- 'test'

# Plot the total number of training observations for each AA in each pos
g <- ggplot(aaCounts.train, aes(x = aa, y=count)) + 
  geom_bar(color = 'black', stat = 'identity') + 
  facet_wrap(~position, scales = 'free_y') + 
  labs(x = 'Amino acid') +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'aa_perPos_training.pdf'),
       height = 6, width = 8)

# Plot the total number of testing observations for each AA in each pos
g <- ggplot(aaCounts.test, aes(x = aa, y=count)) + 
  geom_bar(color = 'black', stat = 'identity') + 
  facet_wrap(~position, scales = 'free_y') + 
  labs(x = 'Amino acid') +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'aa_perPos_testing.pdf'),
       height = 6, width = 8)

# Plot the total number of observations for each AA in each pos (training + testing)
g <- ggplot(rbind(aaCounts.test,aaCounts.train), aes(x = aa, y=count, fill = grp)) + 
  geom_bar(color = 'black', stat = 'identity', position = 'stack') + 
  facet_wrap(~position, scales = 'free_y') + 
  scale_fill_brewer("Group", palette = 'Dark2') + 
  labs(x = 'Amino acid') +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'aa_perPos_both.pdf'),
       height = 6, width = 8)

# Option to remove low IC experimental columns
if (RM_LOW_IC) {
  res <- res[ic.exp >= 0.25]
  RES_DIR <- paste0(RES_DIR, 'rmLowICcols_exp/')
  dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Remove results for proteins where affinity was lost and few 
# k-mers were present to create PWMs
barrera.drop <- unique(s6tab[aff.change == '-' & n.8mers.alt < 10]$prot)
res <- res[!(prot %in% barrera.drop)]
res <- merge(res, matchTab[,c('prot','dset','coreSeq'),with=FALSE])
#res <- res[prot %in% hd1Info$prot]

RES_NONRF <- RES_DIR
if (!INCLUDE_RF_ALI) {
  res <- res[predType != 'rf_preAligned']
} else {
  RES_DIR <- RES_DIR_RF
}
resSumm <- res[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                  numCols = .N),
               by = c('pos','dset','predType')]

# Compare overall per-position performance across the four types of predictions
# on the examples where all models could be used
resSumm.combined <- 
  res[,.(fracAgree = length(which(pcc >= 0.5))/.N,numCols = .N),
      by = c('pos','predType')]
g <- ggplot(resSumm.combined, aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer(palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_combined.pdf'),
       height = 4, width = 6)
g <- ggplot(res, aes(x = pos, y = rmse, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'rmse_boxplots_combined.pdf'),
       height = 4, width = 6)
# Facet per dataset
g <- ggplot(resSumm, aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer(palette = 'Dark2') + 
  facet_wrap(~dset, ncol = length(table(res$dset))) +
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_perDataset.pdf'),
       height = 4, width = 11)
g <- ggplot(res, aes(x = pos, y = rmse, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
  facet_wrap(~dset, ncol = length(table(res$dset))) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'rmse_boxplots_perDataset.pdf'),
       height = 4, width = 11)

# What of we exclude the Chu et al. testing data 
# (I.e., if we're only interested in naturally occurring data)
resSumm.combined.excludeChu2012 <- 
  res[dset != 'Chu2012',.(fracAgree = length(which(pcc >= 0.5))/.N,numCols = .N),
      by = c('pos','predType')]
g <- ggplot(resSumm.combined.excludeChu2012, 
            aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer("Method", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_combined_excludeChu2012.pdf'),
       height = 4, width = 6)
g <- ggplot(res[dset != 'Chu2012'], aes(x = pos, y = rmse, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'rmse_boxplots_combined_excludeChu2012.pdf'),
       height = 4, width = 6)
# Facet per dataset
g <- ggplot(resSumm[dset != 'Chu2012'], aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer("Method", palette = 'Dark2') + 
  facet_wrap(~dset, ncol = length(table(res$dset))) +
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_perDataset_excludeChu2012.pdf'),
       height = 4, width = 11)
g <- ggplot(res[dset != 'Chu2012'], aes(x = pos, y = rmse, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
  facet_wrap(~dset, ncol = length(table(res$dset))) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'rmse_boxplots_perDataset_excludeChu2012.pdf'),
       height = 4, width = 9)

# Compare results for model vs hybrid model+nn
g <- ggplot(resSumm[predType %in% c('model','model+nn')], 
            aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer("Method", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_perDset_comparePredVsTransfer.pdf'),
       height = 4, width = 11)

# Disect the NN vs model predictions for the model+NN groups
disect <- res[predType=='model+nn']
disect$colType <- ifelse(disect$nnbrs == 0, 'model','NN')
disectSumm <- disect[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                        numCols = .N),
                     by = c('pos','dset','colType')]
disectSumm.combined <- disect[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                 numCols = .N),
                              by = c('pos','colType')]
g <- ggplot(disectSumm, aes(x = pos, y = fracAgree, fill = colType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("Column type", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_perDataset.pdf'),
       height = 4, width = 11)
g <- ggplot(disectSumm.combined, aes(x = pos, y = fracAgree, fill = colType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("Column type", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN.pdf'),
       height = 4, width = 6)

# How many do we predict correctly when specificity change occurred in 
# the Barrera dataset?

#### Note ... come back to this
barrera.specChange <- unique(s6tab[spec.change == 'Yes']$prot)

# Plot the model coefficients and get an idea of where things 
# might be going wrong for the chu data
coefDir <- paste0(RES_DIR, 'coefPlots/')
dir.create(coefDir, recursive = TRUE, showWarnings = FALSE)
# How many observations for each aa?
mcoefs$count.train <- sapply(1:nrow(mcoefs), function(i) {
  pos <- mcoefs$aapos[i]; a <- mcoefs$aa[i]
  x <- aaCounts.train[aa == a & position == pos]$count
  if (length(x) == 0) 0 else x
})
wrongPreds <-  
  merge(melt(matchTab.test, id.vars = c('prot','grp'), measure.vars = CORE_POS,
             variable.name = 'position', value.name = 'aa'),
        res[predType == 'model' & pcc < 0.5], by = 'prot')
wrongPreds <- merge(wrongPreds, aaCounts.train, by = c('position','aa'))
wrongPreds$edge <- paste0(wrongPreds$pos,'-',wrongPreds$position)
edges <- sort(unique(paste0(mcoefs$bpos,'-',mcoefs$aapos)))
wrongPreds <- wrongPreds[edge %in% edges]
wrongPreds.count <- wrongPreds[,.(count = .N), by = c('pos','position','aa')]
mcoefs$count.predWrong <- sapply(1:nrow(mcoefs), function(i) {
  ap <- mcoefs$aapos[i]; a <- mcoefs$aa[i]; bp <- mcoefs$bpos[i]
  x <- wrongPreds.count[aa == a & position == ap & pos == bp]$count
  if (length(x) == 0) 0 else x
})
for (bp in 1:6) {
  nr = length(levels(droplevels(factor(mcoefs[bpos == bp]$aapos))))
  g <- ggplot(mcoefs[bpos == bp], aes(x = aa, y = base, fill = coef)) +
    geom_tile(color = 'black') +
    scale_fill_gradient2("Coefficient", low = 'blue', high = 'red')+
    facet_wrap(~aapos, nrow = nr)+
    geom_text(aes(label = count.train), size = 2) +
    theme_classic()
  ggsave(plot = g, file = paste0(coefDir, 'b',bp, '_coefs_nObsAA.pdf'),
         height = 1.2*nr +.5, width = 6)
  
  g <- ggplot(mcoefs[bpos == bp], aes(x = aa, y = base, fill = coef)) +
    geom_tile(color = 'black') +
    scale_fill_gradient2("Coefficient", low = 'blue', high = 'red')+
    facet_wrap(~aapos, nrow = nr)+
    geom_text(aes(label = count.predWrong), size = 2) +
    theme_classic()
  ggsave(plot = g, file = paste0(coefDir, 'b',bp, '_coefs_predWrongCount.pdf'),
         height = 1.2*nr +.5, width = 6)
}

# Info for removing proteins used for training rf_extant or rf_joint
motifInfo.cisbp <- fread('~/research/cisBP/hboxOnlyInfo.txt')
trPubs.rfExtant <- c('Berger08','Zhu09','FlyFactorSurvey')
trSet.rfExtant <- unique(motifInfo.cisbp[MSource_Identifier %in% trPubs.rfExtant]$TF_Name)
#trSet.rfExtant.cores <- unique(res[predType != 'rf_joint' & !(prot %in% trSet.rfExtant))
trSet.rfJoint <- 
  union(unique(motifInfo.cisbp[MSource_Identifier %in% trPubs.rfExtant]$TF_Name),
        unique(matchTab[dset == 'Chu2012']$prot))


# For the Barrera 2016 mutants that vary in only a single base-contacting position,
# how well do we do by transferring columns from the wild-type if they don't
# the column doesn't contact the varied amino acid position, then predicting 
# only the remaining columns?  Do we do better than a de novo prediction from
# the Stormo group?

if (!INCLUDE_RF_ALI) {
  transfer <- fread(paste0(RES_DIR,'pccTable_BarreraMuts_transfer-modelNN.txt'))
  transfer$predType <- 'transfer-model+NN'
  setkey(transfer, prot,pos)
  xferOrder <- transfer$transfer; varPosOrder <- transfer$varPos;
  tmp <- fread(paste0(RES_DIR,'pccTable_BarreraMuts_transfer-model.txt'))
  tmp$predType <- 'transfer-model'
  transfer <- rbind(transfer, tmp)
  tmp <- fread(paste0(stormoFiles[1]))[prot %in% transfer$prot]
  setkey(tmp, prot,pos)
  tmp$transfer <- xferOrder; tmp$varPos <- varPosOrder; tmp$predType <- stormoLabs[1]
  transfer <- rbind(transfer, tmp)
  tmp <- fread(paste0(stormoFiles[2]))[prot %in% transfer$prot]
  setkey(tmp, prot,pos)
  tmp$transfer <- xferOrder; tmp$varPos <- varPosOrder; tmp$predType <- stormoLabs[2]
  transfer <- rbind(transfer, tmp)
  transfer$predType <- 
    factor(transfer$predType, levels = c('transfer-model','transfer-model+NN',stormoLabs))
  transfer$pos <- factor(transfer$pos)
  
  # Toss columns where binding was lost
  transfer <- transfer[!(prot %in% barrera.drop)]
  
  # Get aggregate agreement stats
  transfer$transfer <- ifelse(transfer$transfer, 'transfer possible',
                              'requires predicition')
  transferSumm <- transfer[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                              numCols = .N),
                           by = c('pos','predType')]
  transferSumm.byColType <- transfer[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                        numCols = .N),
                                     by = c('predType','transfer')]
  
  # Compare overall per-position performance
  g <- ggplot(transferSumm, aes(x = pos, y = fracAgree, fill = predType)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_brewer(palette = 'Dark2') + 
    labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_barreraMuts_pcc_fracAgree.pdf'),
         height = 4, width = 6)
  g <- ggplot(transfer, aes(x = pos, y = rmse, fill = predType)) + 
    geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
    #facet_wrap(~dset, ncol = 3) +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_barreraMuts_rmse_boxplots.pdf'),
         height = 4, width = 6)
  
  # Compare per-position performance distinguishing transferred vs predicted columns
  g <- ggplot(transferSumm.byColType, aes(x = transfer, y = fracAgree, fill = predType)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_brewer(palette = 'Dark2') + 
    #facet_wrap(~transfer, ncol = 2) +
    labs(x = "Column type", y = "Columns in agreement (predicted vs. actual)") +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_barreraMuts_pcc_fracAgree_byColType.pdf'),
         height = 4, width = 6)
  g <- ggplot(transfer, aes(x = transfer, y = rmse, fill = predType)) + 
    geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Column type", y = "RMSE (predicted vs. actual)") +
    #facet_wrap(~transfer, ncol = 2) +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_barreraMuts_rmse_boxplots_byColType.pdf'),
         height = 4, width = 6)
  
  
  # What about if we do HD1 transfer for all protein pairs with DBD %ID above a 
  # given threshold?
  hd1Allcmp <- fread(paste0(RES_DIR,'pccTable_allHD1Pairs_transfer-model.txt'))
  hd1Allcmp$predType <- 'transfer_model'
  tmp <- fread(paste0(RES_DIR,'pccTable_allHD1Pairs_transfer-modelNN.txt'))
  tmp$predType <- 'transfer_modelNN'
  hd1Allcmp <- rbind(hd1Allcmp, tmp)
  hd1Allcmp <- merge(hd1Allcmp, matchTab.train[,c('prot','dset')], by.x = 'prot2', by.y = 'prot')
  names(hd1Allcmp)[length(hd1Allcmp)] <- 'dset.prot2'
  hd1Allcmp <- merge(hd1Allcmp, matchTab.test[,c('prot','dset')], by.x = 'prot1', by.y = 'prot')
  names(hd1Allcmp)[length(hd1Allcmp)] <- 'dset.prot1'
  hd1Allcmp$pos <- factor(hd1Allcmp$pos)
  hd1Allcmp$protPair <- paste(hd1Allcmp$prot1, hd1Allcmp$prot2, sep = '-')
  hd1Allcmp$transfer <- ifelse(hd1Allcmp$transfer, 'transfer possible','requires predicition')
  hd1Allcmp <- hd1Allcmp[!(prot1 %in% barrera.drop)]
  hd1Allcmp.addStormo <- hd1Allcmp
  for (i in 1:length(stormoLabs)) {
    tmp <- merge(hd1Allcmp[predType == 'transfer_model',
                           c('prot1','prot2','pos','varPos','transfer','aa',
                              'dset.prot1','dset.prot2')],
          res[predType == stormoLabs[i],c('prot','pos','pcc','rmse','ic.exp',
                                          'ic.pred','predType')],
          by.x = c('prot1','pos'), by.y = c('prot','pos'))
    tmp$protPair <- paste(tmp$prot1, tmp$prot2, sep = '-')
    tmp <- tmp[,names(hd1Allcmp),with = FALSE]
    
    hd1Allcmp.addStormo <- rbind(hd1Allcmp.addStormo[protPair %in% tmp$protPair], tmp)
  }
  # Get aggregate agreement stats
  hd1Allcmp.addStormoSumm <- hd1Allcmp.addStormo[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                              numCols = .N),
                           by = c('pos','predType')]
  hd1Allcmp.addStormoSumm.byColType <- hd1Allcmp.addStormo[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                        numCols = .N),
                                     by = c('predType','transfer')]
  
  # Compare overall per-position performance
  g <- ggplot(hd1Allcmp.addStormoSumm, aes(x = pos, y = fracAgree, fill = predType)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_brewer(palette = 'Dark2') + 
    labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_allHD1_pcc_fracAgree.pdf'),
         height = 4, width = 6)
  g <- ggplot(hd1Allcmp.addStormo, aes(x = pos, y = rmse, fill = predType)) + 
    geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
    #facet_wrap(~dset, ncol = 3) +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_allHD1_rmse_boxplots.pdf'),
         height = 4, width = 6)
  
  # Compare overall per-position performance
  g <- ggplot(hd1Allcmp.addStormoSumm, aes(x = pos, y = fracAgree, fill = predType)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_brewer(palette = 'Dark2') + 
    labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_allHD1_pcc_fracAgree.pdf'),
         height = 4, width = 6)
  g <- ggplot(hd1Allcmp.addStormoSumm[!(predType %in% c('rf_extant','rf_joint'))], 
              aes(x = pos, y = fracAgree, fill = predType)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_brewer(palette = 'Dark2') + 
    labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_allHD1_pcc_fracAgree_noStormo.pdf'),
         height = 4, width = 6)
  g <- ggplot(hd1Allcmp.addStormo, aes(x = pos, y = rmse, fill = predType)) + 
    geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
    #facet_wrap(~dset, ncol = 3) +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_allHD1_rmse_boxplots.pdf'),
         height = 4, width = 6)
  g <- ggplot(hd1Allcmp.addStormo[!(predType %in% c('rf_extant','rf_joint'))], 
              aes(x = pos, y = rmse, fill = predType)) + 
    geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
    #facet_wrap(~dset, ncol = 3) +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_allHD1_rmse_boxplots_noStormo.pdf'),
         height = 4, width = 6)
  
  # Compare per-position performance distinguishing transferred vs predicted columns
  g <- ggplot(hd1Allcmp.addStormoSumm.byColType, aes(x = transfer, y = fracAgree, fill = predType)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_brewer(palette = 'Dark2') + 
    #facet_wrap(~transfer, ncol = 2) +
    labs(x = "Column type", y = "Columns in agreement (predicted vs. actual)") +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_allHD1_pcc_fracAgree_byColType.pdf'),
         height = 4, width = 6)
  g <- ggplot(hd1Allcmp.addStormoSumm.byColType[!(predType %in% c('rf_extant','rf_joint'))], 
              aes(x = transfer, y = fracAgree, fill = predType)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_brewer(palette = 'Dark2') + 
    #facet_wrap(~transfer, ncol = 2) +
    labs(x = "Column type", y = "Columns in agreement (predicted vs. actual)") +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_allHD1_pcc_fracAgree_byColType_noStormo.pdf'),
         height = 4, width = 6)
  g <- ggplot(hd1Allcmp.addStormo, aes(x = transfer, y = rmse, fill = predType)) + 
    geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Column type", y = "RMSE (predicted vs. actual)") +
    #facet_wrap(~transfer, ncol = 2) +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_allHD1_rmse_boxplots_byColType.pdf'),
         height = 4, width = 6)
  g <- ggplot(hd1Allcmp.addStormo[!(predType %in% c('rf_extant','rf_joint'))], 
              aes(x = transfer, y = rmse, fill = predType)) + 
    geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Column type", y = "RMSE (predicted vs. actual)") +
    #facet_wrap(~transfer, ncol = 2) +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'transfer_allHD1_rmse_boxplots_byColType_noStormo.pdf'),
         height = 4, width = 6)
  
  
  # Compare once removing training sets for RF-extant and RF-joint independently
  hd1Allcmp.addStormo.rfExtantOnly <- 
    hd1Allcmp.addStormo[!(prot1 %in% trSet.rfExtant) & predType != 'rf_joint']
  hd1Allcmp.addStormo.rfJointOnly <- 
    hd1Allcmp.addStormo[!(prot1 %in% trSet.rfJoint) & predType != 'rf_extant']  
  # Get aggregate agreement stats
  hd1Allcmp.addStormo.rfExtantOnlySumm <- 
    hd1Allcmp.addStormo.rfExtantOnly[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                                    numCols = .N),
                                                 by = c('pos','predType')]
  hd1Allcmp.addStormo.rfExtantOnly.byColType <-
    hd1Allcmp.addStormo.rfExtantOnly[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                                              numCols = .N),
                                                           by = c('predType','transfer')]
  hd1Allcmp.addStormo.rfJointOnlySumm <- 
    hd1Allcmp.addStormo.rfJointOnly[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                        numCols = .N),
                                     by = c('pos','predType')]
  hd1Allcmp.addStormo.rfJointOnly.byColType <-
    hd1Allcmp.addStormo.rfJointOnly[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                        numCols = .N),
                                     by = c('predType','transfer')]
  
  ### This doesn't look right ...
  
}


hd1Prots <- unique(fread(paste0(RES_NONRF,'pccTable_allHD1Pairs_transfer-modelNN.txt'))$prot1)
# Compare de novo predictions only against the rf_extant predictor, but removing all proteins 
# that it was trained on from the test set
res.noOlap.rfExtant <- res[predType != 'rf_joint' & !(prot %in% trSet.rfExtant)]
hasNeighbors <- unique(res.noOlap.rfExtant[nnbrs > 0]$prot)
res.noOlap.rfExtant.hd1 <- res[predType != 'rf_joint' & !(prot %in% trSet.rfExtant)
                               & prot %in% hd1Prots & prot %in% hasNeighbors]
res.noOlap.rfExtantSumm <- res.noOlap.rfExtant[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                                  numCols = .N),
                                               by = c('pos','predType')]
res.noOlap.rfExtant.hd1Summ <- res.noOlap.rfExtant.hd1[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                                          numCols = .N),
                                                       by = c('pos','predType')]
g <- ggplot(res.noOlap.rfExtantSumm, aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("Column type", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_rfExtantOnly.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfExtantSumm[predType != 'model+nn'], aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic() + 
  theme(legend.position = 'top')
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_rfExtantOnly_excludeNN.pdf'),
       height = 4, width = 4)
g <- ggplot(res.noOlap.rfExtant.hd1Summ, aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("Column type", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_rfExtantOnly_hd1.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfExtant, aes(x = pos, y = pcc, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfExtantOnly.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfExtant[predType != 'model+nn'], aes(x = pos, y = pcc, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("",palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic() + 
  theme(legend.position = 'top')
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfExtantOnly_excludeNN.pdf'),
       height = 4, width = 4)

## Supplemental figure S2, right
s2Tab <- droplevels(res.noOlap.rfExtant[predType != 'model+nn'])
s2Tab$predType <- plyr::mapvalues(s2Tab$predType,from = c('model', 'rf_extant'),to = c('rCLAMPS','rf_extant'))
g <- ggplot(s2Tab, aes(x = pos, y = pcc, grp = predType)) + 
  geom_lv(color = 'black', outlier.size = 1, aes(fill = predType), alpha = .7) + 
  scale_fill_brewer("", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  scale_y_continuous(limits = c(-1,1)) +
  theme_classic() + 
  theme(legend.position = 'top')
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfExtantOnly_excludeNN_boxen.pdf'),
       height = 4, width = 4)

g <- ggplot(res.noOlap.rfExtant.hd1, aes(x = pos, y = pcc, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfExtantOnly_hd1.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfExtant, aes(x = pos, y = rmse, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  #geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'rmse_boxplots_rfExtantOnly.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfExtant.hd1, aes(x = pos, y = rmse, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  #geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'rmse_boxplots_rfExtantOnly_hd1.pdf'),
       height = 4, width = 6)

# Compare de novo predictions only against the rf_extant predictor, but removing all proteins 
# that it was trained on from the test set
res.noOlap.rfJoint <- res[predType != 'rf_extant' & !(prot %in% trSet.rfJoint)]
hasNeighbors <- unique(res.noOlap.rfJoint[nnbrs > 0]$prot)
res.noOlap.rfJoint.hd1 <- res[predType != 'rf_extant' & !(prot %in% trSet.rfJoint) 
                              & prot %in% hd1Prots & prot %in% hasNeighbors]
res.noOlap.rfJointSumm <- res.noOlap.rfJoint[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                                  numCols = .N),
                                               by = c('pos','predType')]
res.noOlap.rfJoint.hd1Summ <- res.noOlap.rfJoint.hd1[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                                        numCols = .N),
                                                     by = c('pos','predType')]

g <- ggplot(res.noOlap.rfJointSumm, aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("Column type", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_rfJointOnly.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfJointSumm[predType != 'model+nn'], aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic() + 
  theme(legend.position = 'top')
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_rfJointOnly_excludeNN.pdf'),
       height = 4, width = 4)

## S3 left
s3left <- droplevels(res.noOlap.rfJointSumm[predType != 'model+nn'])
s3left$predType <- plyr::mapvalues(s3left$predType,from = c('model', 'rf_joint'),to = c('rCLAMPS','rf_joint'))
g <- ggplot(s3left, aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic() + 
  theme(legend.position = 'top')
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_rfJointOnly_excludeNN_rCLAMPS.pdf'),
       height = 4, width = 4)

## S2 left
s2left <- droplevels(res.noOlap.rfExtantSumm[predType != 'model+nn'])
s2left$predType <- plyr::mapvalues(s2left$predType,from = c('model', 'rf_extant'),to = c('rCLAMPS','rf_extant'))
g <- ggplot(s2left, aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic() + 
  theme(legend.position = 'top')
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_rfExtantOnly_excludeNN_rCLAMPS.pdf'),
       height = 4, width = 4)

g <- ggplot(res.noOlap.rfJoint.hd1Summ, aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("Column type", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_rfJointOnly_hd1.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfJoint, aes(x = pos, y = pcc, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfJointOnly.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfJoint[predType != 'model+nn'], aes(x = pos, y = pcc, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic() +
  theme(legend.position = 'top')
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfJointOnly_excludeNN.pdf'),
       height = 4, width = 4)

## Supplemental figure S3, right
s3Tab <- droplevels(res.noOlap.rfJoint[predType != 'model+nn'])
s3Tab$predType <- plyr::mapvalues(s3Tab$predType,from = c('model', 'rf_joint'),to = c('rCLAMPS','rf_joint'))
g <- ggplot(s3Tab, aes(x = pos, y = pcc, grp = predType)) + 
  geom_lv(color = 'black', outlier.size = 1, aes(fill = predType), alpha = .7) + 
  scale_fill_brewer("", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  scale_y_continuous(limits = c(-1,1)) +
  theme_classic() + 
  theme(legend.position = 'top')
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfJointOnly_excludeNN_boxen.pdf'),
       height = 4, width = 4)

g <- ggplot(res.noOlap.rfJoint.hd1, aes(x = pos, y = pcc, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC between predicted and actual") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfJointOnly_hd1.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfJoint, aes(x = pos, y = rmse, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  #geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'rmse_boxplots_rfJointOnly.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfJoint.hd1, aes(x = pos, y = rmse, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  #geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'rmse_boxplots_rfJointOnly_hd1.pdf'),
       height = 4, width = 6)

# For the fair comparisons, how do the predicitons break down in 
# terms of NN (transferred) vs. model-predicted columns
res.noOlap.rfExtant.hd1 <- merge(res.noOlap.rfExtant.hd1,
                                 res.noOlap.rfExtant.hd1[predType == 'model+nn',
                                                         .(transfer = nnbrs > 0),by=c('prot','pos')],
                                 by = c('prot','pos'))
res.noOlap.rfExtant.hd1$transfer <- ifelse(res.noOlap.rfExtant.hd1$transfer, "transfer\npossible", "requires\npredicition")
res.noOlap.rfExtant.hd1Summ.byColType <- res.noOlap.rfExtant.hd1[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                                                    numCols = .N),
                                                                 by = c('predType','transfer')]
res.noOlap.rfJoint.hd1 <- merge(res.noOlap.rfJoint.hd1,
                                res.noOlap.rfJoint.hd1[predType == 'model+nn',
                                                       .(transfer = nnbrs > 0),by=c('prot','pos')],
                                by = c('prot','pos'))
res.noOlap.rfJoint.hd1$transfer <- ifelse(res.noOlap.rfJoint.hd1$transfer, "transfer\npossible", "requires\npredicition")
res.noOlap.rfJoint.hd1Summ.byColType <- res.noOlap.rfJoint.hd1[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                                                  numCols = .N),
                                                               by = c('predType','transfer')]

# Compare per-position performance distinguishing transferred vs predicted columns
g <- ggplot(res.noOlap.rfExtant.hd1Summ.byColType, aes(x = transfer, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer(palette = 'Dark2') + 
  #facet_wrap(~transfer, ncol = 2) +
  labs(x = "Column type", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_rfExtantOnly_hd1_byColType.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfExtant.hd1, aes(x = transfer, y = pcc, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Column type", y = "PCC between predicted and actual") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfExtantOnly_hd1_byColType.pdf'),
       height = 4, width = 4)

#### For Figure 5
library(dplyr)
tmp <- melt(res.noOlap.rfExtant.hd1[transfer == 'transfer\npossible'],
            id.vars = c('prot', 'pos', 'predType', 'dset','coreSeq'), measure.vars = c('pcc','rmse'),
            value.name = 'score', variable.name = 'measure')
tmp$predType <- factor(plyr::mapvalues(tmp$predType, from = c('model','model+nn','rf_extant'),
                                to = c('model','transfer','rf_extant')),
                       levels = c('transfer', 'model','rf_extant'))
tmp$measure <- factor(tmp$measure, levels = c('rmse','pcc'))

g <- ggplot(tmp, aes(x = predType, y = score)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  #geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Method", y = "Agreement (predicted vs. actual)") +
  facet_wrap(~measure, ncol = 2, scales = 'free_y') +
  theme_bw() + 
  theme(axis.title.x = element_blank())
ggsave(plot = g, file = paste0(RES_DIR, 'boxplots_rfExtantOnly_hd1_pcc-rmse_combined.pdf'),
       height = 4, width = 8)

g <- ggplot(res.noOlap.rfExtant.hd1, aes(x = transfer, y = rmse, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  #geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Column type", y = "RMSE (predicted vs. actual)") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'rmse_boxplots_rfExtantOnly_hd1_byColType.pdf'),
       height = 4, width = 4)


g <- ggplot(res.noOlap.rfJoint.hd1Summ.byColType, aes(x = transfer, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer(palette = 'Dark2') + 
  #facet_wrap(~transfer, ncol = 2) +
  labs(x = "Column type", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_compareModelVsNN_rfJointOnly_hd1_byColType.pdf'),
       height = 4, width = 6)
g <- ggplot(res.noOlap.rfJoint.hd1, aes(x = transfer, y = pcc, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Column type", y = "PCC between predicted and actual") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfJointOnly_hd1_byColType.pdf'),
       height = 4, width = 6)

res.noOlap.rfJoint.hd1$predType <- plyr::mapvalues(res.noOlap.rfJoint.hd1$predType,
                                                   from = c('model', 'model+nn','rf_extant','rf_joint'),
                                                   to = c('de novo', 'hybrid','rf_extant','rf_joint'))
g <- ggplot(res.noOlap.rfJoint.hd1[transfer == 'transfer\npossible' & predType %in% c('hybrid','rf_joint')], 
            aes(x = predType, y = pcc, fill = predType)) + 
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  #geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  #scale_fill_brewer("Method", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Method", y = "PCC between predicted and actual") +
  #facet_wrap(~dset, ncol = 3) +
  scale_y_continuous(limits = c(-1,1)) +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfJointOnly_hd1_byColType_forSupplement.pdf'),
       height = 5, width = 5)

## Figure 5 - plot the overall comparison of transfer vs OUR model
res.fig5 <- merge(res,
                  res[predType == 'model+nn',.(transfer = nnbrs > 0),by=c('prot','pos')],
                  by = c('prot','pos'))
res.fig5$transfer <- ifelse(res.fig5$transfer, "transfer", "de novo")
nearMuts <- unique(res.fig5[transfer == 'transfer']$prot)
res.fig5 <- res.fig5[prot %in% nearMuts]
res.fig5$agree <- ifelse(res.fig5$pcc >= 0.5, TRUE, FALSE)
res.fig5$predType <- plyr::mapvalues(res.fig5$predType,
                                     from = c('model', 'model+nn','rf_extant','rf_joint'),
                                     to = c('model', 'hybrid','rf_extant','rf_joint'))
res.fig5$predType <- factor(res.fig5$predType,
                            levels = c('hybrid','model','rf_extant','rf_joint'))


# Figure 5 boxen plots for hybrid approach only
g <- ggplot(res.fig5[predType == 'hybrid'], aes(x = transfer, y = pcc)) + 
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  scale_y_continuous(limits = c(-1,1)) +
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Column type", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxen_plots_hybridApproach_byColType.pdf'),
       height = 4, width = 3)

g <- ggplot(res.fig5[predType == 'hybrid'], aes(x = transfer, y = pcc)) + 
  geom_boxplot(color = 'black', outlier.size = 0.6) +
  scale_y_continuous(limits = c(-1,1)) +
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Column type", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_hybridApproach_byColType.pdf'),
       height = 4, width = 3)

# Repeat with RMSE and pcc side by side
tmp <- melt(res.fig5[predType == 'hybrid'], 
            id.vars = c('prot','pos','predType','transfer'),
            measure.vars = c('rmse','pcc')) 
g <- ggplot(tmp, aes(x = transfer, y = value)) + 
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  #geom_hline(data = tmp[variable == 'pcc'], yintercept = 0.5, lty = 'dashed') +
  labs(x = "Column type", y = "Agreement (predicted vs. actual)") +
  facet_wrap(~variable, scales = 'free_y') + 
  theme_bw()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxen_plots_hybridApproach_byColType_addRMSE.pdf'),
       height = 4, width = 4)

# Figure 5 boxen plots rf_extant comparison
g <- ggplot(res.fig5[!(prot %in% trSet.rfExtant) &
                       predType %in% c('model','hybrid','rf_extant') & transfer == 'transfer'], 
            aes(x = predType, y = pcc)) + 
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  scale_y_continuous(limits = c(-1,1)) +
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Method", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxen_plots_rfExtantHD1_transferOnly.pdf'),
       height = 4, width = 3)

g <- ggplot(res.fig5[!(prot %in% trSet.rfExtant) &
                       predType %in% c('model','hybrid','rf_extant') & transfer == 'transfer'], 
            aes(x = predType, y = pcc)) + 
  geom_boxplot(color = 'black', outlier.size = 0.6) + 
  scale_y_continuous(limits = c(-1,1)) +
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Method", y = "PCC between predicted and actual") +
  theme_classic()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxplots_rfExtantHD1_transferOnly.pdf'),
       height = 4, width = 3)

# Repeat with RMSE and pcc side by side
tmp <- melt(res.fig5[!(prot %in% trSet.rfExtant) &
                       predType %in% c('model','hybrid','rf_extant') & transfer == 'transfer'], 
            id.vars = c('prot','pos','predType','transfer'),
            measure.vars = c('rmse','pcc')) 
g <- ggplot(tmp, aes(x = predType, y = value)) + 
  geom_lv(color = 'black', fill = "gray20", outlier.size = 1, alpha = 0.3) + 
  #geom_hline(data = tmp[variable == 'pcc'], yintercept = 0.5, lty = 'dashed') +
  labs(x = "Method", y = "Agreement (predicted vs. actual)") +
  facet_wrap(~variable, scales = 'free_y') + 
  theme_bw()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_boxen_plots_rfExtantHD1_transferOnly_addRMSE.pdf'),
       height = 4, width = 5)

# Number of columns transferred vs. number of columns predicted de novo by model+nn
table(res.fig5[predType == 'model+nn']$transfer)
# Fraction of model+nn columns predicted correctly when trasferred vs. predicted de novo
table(res.fig5[predType == 'model+nn']$transfer, res.fig5[predType == 'model+nn']$agree)


print(paste("# Proteins in rf_extant comparison:",length(unique(res.noOlap.rfExtant$prot))))
print(paste("# Core seqs in rf_extant comparison:",length(unique(res.noOlap.rfExtant$coreSeq))))
print(paste("# Proteins in rf_extant HD1 comparison:",length(unique(res.noOlap.rfExtant.hd1$prot))))
print(paste("# Core seqs in rf_extant HD1 comparison:",length(unique(res.noOlap.rfExtant.hd1$coreSeq))))
print(paste("# Proteins in rf_joint comparison:",length(unique(res.noOlap.rfJoint$prot))))
print(paste("# Core seqs in rf_joint comparison:",length(unique(res.noOlap.rfJoint$coreSeq))))
print(paste("# Proteins in rf_joint comparison:",length(unique(res.noOlap.rfJoint.hd1$prot))))
print(paste("# Core seqs in rf_joint comparison:",length(unique(res.noOlap.rfJoint.hd1$coreSeq))))

# Save some of the outputs from this analysis
write.table(res.noOlap.rfJoint.hd1,file = paste0(RES_DIR,'table_res.noOlap.rfJoint.hd1.txt'),
            quote = FALSE, sep = '\t', row.names = FALSE)
write.table(res.noOlap.rfJoint,file = paste0(RES_DIR,'table_res.noOlap.rfJoint.txt'),
            quote = FALSE, sep = '\t', row.names = FALSE)
write.table(res.noOlap.rfExtant.hd1,file = paste0(RES_DIR,'table_res.noOlap.rfExtant.hd1.txt'),
            quote = FALSE, sep = '\t', row.names = FALSE)
write.table(res.noOlap.rfExtant,file = paste0(RES_DIR,'table_res.noOlap.rfExtant.txt'),
            quote = FALSE, sep = '\t', row.names = FALSE)
write.table(res,file = paste0(RES_DIR,'table_res.txt'),quote = FALSE, 
            sep = '\t', row.names = FALSE)

# Find the fraction with PCC >= x for the two main test datasets
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
  g <- ggplot(fracPCCge.list[[dset.name[i]]], aes(x = pcc, y = fracPCCge)) + 
    geom_line(size = 0.75, aes(col = colType)) + 
    geom_vline(xintercept = 0.5, lty = 'dashed') +
    scale_color_brewer("", palette = 'Dark2') + 
    labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
    theme_bw()+ 
    theme(legend.position = 'top')
  #ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracGE_',dset.name[i],'.pdf'),height = 4, width = 6)
  ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracGE_',dset.name[i],'.pdf'),height = 4.5, width = 4)
  
  # Exclude model+nn
  g <- ggplot(fracPCCge.list[[dset.name[i]]][colType != 'model+nn'], aes(x = pcc, y = fracPCCge)) + 
    geom_line(size = 0.75, aes(col = colType)) + 
    geom_vline(xintercept = 0.5, lty = 'dashed') +
    scale_color_brewer("", palette = 'Dark2') + 
    labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
    theme_bw()+ 
    theme(legend.position = 'top')
  #ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracGE_',dset.name[i],'.pdf'),height = 4, width = 6)
  ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracGE_',dset.name[i],'_encludeNN.pdf'),height = 4.5, width = 4)
}

# Plot faceted version
rfJointPCCge <- fracPCCge.list[['rf_joint']]
rfJointPCCge$dset <- 'rf_joint'
rfExtantPCCge <- fracPCCge.list[['rf_extant']]
rfExtantPCCge$dset <- 'rf_extant'
pccGE.tab <- rbind(rfExtantPCCge, rfJointPCCge)
#pccGE.tab$dset <- factor(pccGE.tab$dset, levels = c('rf_joint','rf_extant','model','model+nn'))
g <- ggplot(droplevels(pccGE.tab), aes(x = pcc, y = fracPCCge)) + 
  geom_line(size = 0.75, aes(col = colType)) + 
  geom_vline(xintercept = 0.5, lty = 'dashed') +
  scale_color_brewer("", palette = 'Dark2') + 
  labs(x = "PCC", y = "Fraction of columns with PCC >= x") +
  facet_wrap(~dset, scales = 'free_y')+
  theme_bw()
ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracGE_combined_faceted.pdf'),height = 4, width = 10)

if (FALSE) {
  # Exclude columns where transfer prediction could not be made
  resSumm.nnbrsGT0 <- res[nnbrs > 0,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                      numCols = .N),
                          by = c('pos','dset','predType')]
  
  # How many columns have at least one neighbor for transfer (in each base position)
  g <- ggplot(resSumm.nnbrsGT0[predType == 'transfer'],
              aes(x = pos, y = numCols, fill = dset)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
    scale_fill_brewer(palette = 'Dark2') + 
    labs(x = "Binding site position", y = "Number of columns with >= 1 neighbor") +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'numColumns_with_nnbrsGT0.pdf'),
         height = 4, width = 6)
  
  # Make a plot for transfer method, excluding the columns where prediction
  # was required due to lack of neighbors
  g <- ggplot(resSumm.nnbrsGT0[predType == 'transfer'],
              aes(x = pos, y = fracAgree, fill = dset)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
    scale_y_continuous(limits = c(0,1)) +
    geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
    scale_fill_brewer("Dataset", palette = 'Dark2') + 
    labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_perDset_transfer_nnbrsGT0.pdf'),
         height = 4, width = 6)
  
  # Compare results for transfer method vs prediction, excluding the columns where prediction
  # was required due to lack of neighbors
  g <- ggplot(resSumm.nnbrsGT0, aes(x = pos, y = fracAgree, fill = predType)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
    scale_y_continuous(limits = c(0,1)) +
    geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
    facet_wrap(~dset, ncol = 3) +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_perDset_comparePredVsTransfer_nnbrsGT0.pdf'),
         height = 4, width = 11)
  
  res$nnbrs.predType <- ifelse(res$nnbrs == 0, 'prediction','transfer')
  # Exclude columns where transfer prediction could not be made
  resSumm.nnbrs.predType <- res[predType=='transfer',
                                .(fracAgree = length(which(pcc >= 0.5))/.N,
                                  numCols = .N),
                                by = c('pos','dset','nnbrs.predType')]
  g <- ggplot(resSumm.nnbrs.predType, aes(x = pos, y = fracAgree, fill = nnbrs.predType)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
    scale_y_continuous(limits = c(0,1)) +
    geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
    facet_wrap(~dset, ncol = 3) +
    theme_classic()
  ggsave(plot = g, file = paste0(RES_DIR, 'pcc_fracAgree_comparePredVsTransfer.pdf'),
         height = 4, width = 11)
}

