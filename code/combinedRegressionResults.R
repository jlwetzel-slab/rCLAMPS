library(data.table)
library(ggplot2)
library(RColorBrewer)
rm(list = ls())

TAB_DIRS <- c('randomForest_regression_scaled50/allEdges/',
              'randomForest_regression_scaled50/structEdges/',
              'randomForest_regression_scaled50_optuna/',
              'gradientBoosting_regression_scaled50_optuna/')
TAB_DIRS <- paste0('../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/',
                   TAB_DIRS)
TAB_DIRS_TRAINFIT <- TAB_DIRS
TAB_DIRS <- paste0(TAB_DIRS, 'testFit/')
MODEL_LABS <- c('rf_allEdges','rf_structEdges','rf_optuna_best','gb_optuna_best')
TAB_FILES <- c('res','res.noOlap.rfExtant','res.noOlap.rfJoint',
               'res.noOlap.rfExtant.hd1','res.noOlap.rfJoint.hd1')
TAB_FILES_TRAINFIT <- rep('pccTable_underS.txt', length(TAB_DIRS))
TAB_LABS <- c('combined','rfExtantOnly','rfJointOnly','rfExtantOnly_hd1','rfJointOnly_hd1')
IC_THRESH <- 0.25
PLOT_DIR <- '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/'
if (IC_THRESH > 0) {
  PLOT_DIR <- paste0(PLOT_DIR, 'regressionResults_combined_ICge',IC_THRESH,'/')
} else {
  PLOT_DIR <- paste0(PLOT_DIR, 'regressionResults_combined/')
}
dir.create(PLOT_DIR,showWarnings = FALSE,recursive = TRUE)

# Make plots for the in sample fit (i.e., training data fit)
tabs_train <- fread('../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/pccTable_underS.txt')
tabs_train$predType <- 'glm'
for (j in 1:length(TAB_DIRS)) {
  modelLab <- MODEL_LABS[j]
  tabDir <- TAB_DIRS_TRAINFIT[j]
  tmp <- fread(paste0(tabDir,TAB_FILES_TRAINFIT[j]))
  tmp$predType <- modelLab
  tabs_train <- rbind(tabs_train, tmp)
}
tabs_train$pos <- factor(tabs_train$pos)
tabs_train$predType <- factor(tabs_train$predType,
                              levels = c('glm',MODEL_LABS))
fracCorrect <- tabs_train[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                             numCols = .N),by = c(c('pos','predType'))]

g <- ggplot(fracCorrect, aes(x = pos, y = fracAgree, fill = predType)) + 
  geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
  scale_y_continuous(limits = c(0,1)) +
  #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
  scale_fill_brewer("Method", palette = 'Dark2') + 
  labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
  theme_classic()
ggsave(plot = g, file = paste0(PLOT_DIR, 'trainingFit_pcc_fracAgree.pdf'),
       height = 4, width = 6)
g <- ggplot(tabs_train, aes(x = pos, y = pcc, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "PCC (predicted vs. actual)") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(PLOT_DIR,'trainingFit_pcc_boxplots.pdf'),
       height = 4, width = 6)
g <- ggplot(tabs_train, aes(x = pos, y = rmse, fill = predType)) + 
  geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  #geom_hline(yintercept = 0.5, lty = 'dashed') +
  labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
  #facet_wrap(~dset, ncol = 3) +
  theme_classic()
ggsave(plot = g, file = paste0(PLOT_DIR, 'trainingFit_rmse_boxplots.pdf'),
       height = 4, width = 6)

# Make plots for various comparisons out of sample (testing fits)
tabs <- NULL
for (i in 1:length(TAB_FILES)) {
  tabFile <- TAB_FILES[i]
  tabLab <- TAB_LABS[i]
  print(tabLab)
  thisTab <- NULL
  for (j in 1:length(TAB_DIRS)) {
    modelLab <- MODEL_LABS[j]
    tabDir <- TAB_DIRS[j]
    tmp <- fread(paste0(tabDir,'table_',tabFile,'.txt'))
    tmp$predType <- ifelse(tmp$predType == 'rf_preAligned', modelLab, tmp$predType)
    tmp$predType <- ifelse(tmp$predType == 'model', 'glm', tmp$predType)
    tmp$predType <- ifelse(tmp$predType == 'model+nn', 'glm+nn', tmp$predType)
    if (j == 1) {
      thisTab <- rbind(thisTab, tmp[predType %in% c('glm','glm+nn','rf_extant','rf_joint')])
    }
    thisTab <- rbind(thisTab,tmp[predType == modelLab])
  }
  thisTab$ic.pcc <- 0.5*thisTab$ic.exp*thisTab$pcc
  thisTab$pos <- factor(thisTab$pos)
  thisTab$predType <- factor(thisTab$predType, 
                             levels = c('glm','glm+nn',MODEL_LABS,'rf_extant','rf_joint'))
  thisTab$dset <- factor(thisTab$dset)
  thisTab <- thisTab[ic.exp >= IC_THRESH]
  tabs[[tabLab]] <- thisTab
  print(paste("Test set:", tabLab))
  print(paste("# distinct proteins used in testing:",length(unique(thisTab$prot))))
  print(paste("# Binding distinct residue combinations used in testing:",length(unique(thisTab$coreSeq))))
  fracCorrect <- thisTab[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                            numCols = .N),by = c(c('pos','predType'))]
  g <- ggplot(fracCorrect, aes(x = pos, y = fracAgree, fill = predType)) + 
    geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
    scale_y_continuous(limits = c(0,1)) +
    #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Binding site position", y = "Columns in agreement (predicted vs. actual)") +
    theme_classic()
  ggsave(plot = g, file = paste0(PLOT_DIR, tabLab, '_pcc_fracAgree.pdf'),
         height = 4, width = 6)
  g <- ggplot(thisTab, aes(x = pos, y = pcc, fill = predType)) + 
    geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
    scale_fill_brewer("Method", palette = 'Dark2') + 
    geom_hline(yintercept = 0.5, lty = 'dashed') +
    labs(x = "Binding site position", y = "PCC (predicted vs. actual)") +
    #facet_wrap(~dset, ncol = 3) +
    theme_classic()
  ggsave(plot = g, file = paste0(PLOT_DIR, tabLab, '_pcc_boxplots.pdf'),
         height = 4, width = 6)
  g <- ggplot(thisTab, aes(x = pos, y = ic.pcc, fill = predType)) + 
    geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
    scale_fill_brewer("Method", palette = 'Dark2') + 
    labs(x = "Binding site position", y = "IC-PCC (predicted vs. actual)") +
    #facet_wrap(~dset, ncol = 3) +
    theme_classic()
  ggsave(plot = g, file = paste0(PLOT_DIR, tabLab, '_ic-pcc_boxplots.pdf'),
         height = 4, width = 6)
  g <- ggplot(thisTab, aes(x = pos, y = rmse, fill = predType)) + 
    geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
    scale_fill_brewer("Method", palette = 'Dark2') + 
    #geom_hline(yintercept = 0.5, lty = 'dashed') +
    labs(x = "Binding site position", y = "RMSE (predicted vs. actual)") +
    #facet_wrap(~dset, ncol = 3) +
    theme_classic()
  ggsave(plot = g, file = paste0(PLOT_DIR, tabLab, '_rmse_boxplots.pdf'),
         height = 4, width = 6)
  if (length(grep('hd1', tabLab)) > 0) {
    fracByColType <- thisTab[,.(fracAgree = length(which(pcc >= 0.5))/.N,
                                numCols = .N),
                             by = c('predType','transfer')]
    g <- ggplot(fracByColType, aes(x = transfer, y = fracAgree, fill = predType)) + 
      geom_bar(color = 'black', width = 0.6, stat = 'identity', position = 'dodge') + 
      scale_y_continuous(limits = c(0,1)) +
      #geom_text(aes(label=numCols), position=position_dodge(width=0.6), vjust=-0.25,size=2)+
      scale_fill_brewer("Column type", palette = 'Dark2') + 
      labs(x = "Column type", y = "Columns in agreement (predicted vs. actual)") +
      theme_classic()
    ggsave(plot = g, file = paste0(PLOT_DIR, tabLab, '_pcc_fracAgree_byColType.pdf'),
           height = 4, width = 6)
    g <- ggplot(thisTab, aes(x = transfer, y = pcc, fill = predType)) + 
      geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
      scale_fill_brewer("Method", palette = 'Dark2') + 
      geom_hline(yintercept = 0.5, lty = 'dashed') +
      labs(x = "Column type", y = "PCC (predicted vs. actual)") +
      #facet_wrap(~dset, ncol = 3) +
      theme_classic()
    ggsave(plot = g, file = paste0(PLOT_DIR, tabLab, '_pcc_boxplots_byColType.pdf'),
           height = 4, width = 6)
    g <- ggplot(thisTab, aes(x = transfer, y = ic.pcc, fill = predType)) + 
      geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
      scale_fill_brewer("Method", palette = 'Dark2') + 
      labs(x = "Column type", y = "IC-PCC (predicted vs. actual)") +
      #facet_wrap(~dset, ncol = 3) +
      theme_classic()
    ggsave(plot = g, file = paste0(PLOT_DIR, tabLab, '_ic-pcc_boxplots_byColType.pdf'),
           height = 4, width = 6)
    g <- ggplot(thisTab, aes(x = transfer, y = rmse, fill = predType)) + 
      geom_boxplot(color = 'black', width = 0.6, position = 'dodge', outlier.size = 0.1) + 
      scale_fill_brewer("Method", palette = 'Dark2') + 
      labs(x = "Column type", y = "RMSE (predicted vs. actual)") +
      #facet_wrap(~dset, ncol = 3) +
      theme_classic()
    ggsave(plot = g, file = paste0(PLOT_DIR, tabLab, '_rmse_boxplots_byColType.pdf'),
           height = 4, width = 6)
  }
}

# When we predict a position incorrectly, are there other examples for 
# the same core sequence, and if so do we ever get the corresponding position
# right for the other examples?
predPWMs.glm <- fread('../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/transfer_test/fracID_fullDBD/pwms_glmOnly_testSet.txt')
predPWMs.glm$bpos <- factor(predPWMs.glm$bpos)
coreCount <- tabs$combined[predType == 'glm', .(nCore = .N/6),by = c('coreSeq')]
glmWrong <- merge(tabs$combined[predType == 'glm' & pcc < 0.5], coreCount, by = 'coreSeq')
wrongVsRight <- glmWrong[,.(nCore = nCore[1], nWrong = .N,
                            chuOnly = length(dset) == 1 && dset == 'Chu2012'),
                            by = c('pos','coreSeq')]
g <- ggplot(wrongVsRight, aes(x = nCore, y = nWrong)) + 
  geom_jitter(shape = 21, aes(color = chuOnly)) + 
  geom_abline(intercept = 0, slope = 1, lty = 'dashed') +
  scale_x_continuous(limits = c(-0.5,max(wrongVsRight$nCore))) + 
  scale_y_continuous(limits = c(-0.5,max(wrongVsRight$nCore))) + 
  labs(x = "Multiplicity of residue combination",
       y = "Number of times column predicted incorrectly") + 
  scale_fill_brewer("Method", palette = 'Dark2') + 
  facet_wrap(~pos) +
  theme_classic()
ggsave(plot = g, file = paste0(PLOT_DIR, 'glm_timesWrongPerCore.pdf'),
       height = 4, width = 6)

