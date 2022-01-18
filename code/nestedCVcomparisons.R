library(data.table)
library(ggplot2)
library(RColorBrewer)
rm(list = ls())

getResTabs <- function(dirs, labs, n_cv, trainType) {
  res <- NULL
  for (i in 1:length(dirs)) {
    for (j in 0:(n_cv-1)) {
      tmp <- fread(paste0(dirs[i],'pccTable_test_pred_fold',j,'.txt'))
      tmp$modelType <- factor(labs[i])
      tmp$fold <- factor(j)
      tmp$trainType <- factor(trainType)
      res <- rbind(res, tmp)
    }
  }
  res
}

N_CV <- 10
GLM_UNSCALED <- '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15/'
GLM_SCALED <- '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50/'
GLM_UNSCALED_RF_UNSCALED <- 
  paste0(paste0(GLM_UNSCALED,paste0('randomForest_regression/',
                                    c('allEdges/','structEdges/'))),'nestedCV_fit/')
GLM_UNSCALED_RF_SCALED <- 
  paste0(paste0(GLM_UNSCALED,paste0('randomForest_regression_scaled50/',
                                  c('allEdges/','structEdges/'))),'nestedCV_fit/')
GLM_SCALED_RF_SCALED <- 
  paste0(paste0(GLM_SCALED,paste0('randomForest_regression_scaled50/',
                                  c('allEdges/','structEdges/'))),'nestedCV_fit/')
RF_TYPES <- c('complete core', 'contact-aware')
OUT_DIR <- paste0(GLM_SCALED,'nestedCVcomparisons/')

# Compare nested CV fit for scaled vs unscaled data based on the unscaled GLM results
res <- getResTabs(GLM_UNSCALED_RF_UNSCALED,RF_TYPES,N_CV,'rf_unscaled')
res <- rbind(res,getResTabs(GLM_UNSCALED_RF_SCALED,RF_TYPES,N_CV,'rf_scaled'))
res$pos <- factor(res$pos)

# Box plots
# Per base position
g <- ggplot(res, aes(x = pos, y = pcc, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "PCC") +
  geom_hline(yintercept = 0.5, lty='dashed') +
  facet_wrap(~trainType)+
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccBoxes_glm_unscaled_perPos.pdf'),
       width = 8, height = 3.25)
g <- ggplot(res, aes(x = pos, y = rmse, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "rmse") +
  facet_wrap(~trainType)+
  theme_classic()  
ggsave(plot = g, file = paste0(OUT_DIR, 'rmseBoxes_glm_unscaled_perPos.pdf'),
       width = 8, height = 3.25)
# Per fold
g <- ggplot(res, aes(x = fold, y = pcc, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Fold number", y = "PCC") +
  geom_hline(yintercept = 0.5, lty='dashed') +
  facet_wrap(~trainType)+
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccBoxes_glm_unscaled_perFold.pdf'),
       width = 8, height = 3.25)
g <- ggplot(res, aes(x = fold, y = rmse, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Fold number", y = "rmse") +
  facet_wrap(~trainType)+
  theme_classic()  
ggsave(plot = g, file = paste0(OUT_DIR, 'rmseBoxes_glm_unscaled_perFold.pdf'),
       width = 8, height = 3.25)
  
# Aggregate plots
# Per base position
tmp <- res[,.(fracAgree = length(which(pcc>=0.5))/.N),by = c('pos','modelType','trainType')]
g <- ggplot(tmp, aes(x = pos, y = fracAgree, fill = modelType)) +
  geom_bar(color = 'black', position = 'dodge', stat = 'identity') + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "Fraction of coumns in agreement") +
  facet_wrap(~trainType)+
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccFracAgree_glm_unscaled_perPos.pdf'),
       width = 8, height = 3.25)
tmp <- res[,.(fracAgree = length(which(pcc>=0.5))/.N),by = c('fold','modelType','trainType')]
g <- ggplot(tmp, aes(x = fold, y = fracAgree, fill = modelType)) +
  geom_bar(color = 'black', position = 'dodge', stat = 'identity') + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Fold number", y = "Fraction of coumns in agreement") +
  facet_wrap(~trainType)+
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccFracAgree_glm_unscaled_perFold.pdf'),
       width = 8, height = 3.25)

# Per fold
g <- ggplot(res, aes(x = fold, y = pcc, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "PCC") +
  geom_hline(yintercept = 0.5, lty='dashed') +
  facet_wrap(~trainType)+
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccBoxes_glm_unscaled_perFold.pdf'),
       width = 8, height = 3.25)
g <- ggplot(res, aes(x = fold, y = rmse, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "rmse") +
  facet_wrap(~trainType)+
  theme_classic()  
ggsave(plot = g, file = paste0(OUT_DIR, 'rmseBoxes_glm_unscaled_perFold.pdf'),
       width = 8, height = 3.25)
res.unscaled <- res

# Now look at what happens with the scaled GLM data
res <- getResTabs(GLM_SCALED_RF_SCALED,RF_TYPES,N_CV,'rf_scaled')
res$pos <- factor(res$pos)
# Box plots
# Per base position
g <- ggplot(res, aes(x = pos, y = pcc, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "PCC") +
  geom_hline(yintercept = 0.5, lty='dashed') +
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccBoxes_glm_scaled_perPos.pdf'),
       width = 6, height = 3.25)
g <- ggplot(res, aes(x = pos, y = rmse, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "rmse") +
  theme_classic()  
ggsave(plot = g, file = paste0(OUT_DIR, 'rmseBoxes_glm_scaled_perPos.pdf'),
       width = 6, height = 3.25)
# Per fold
g <- ggplot(res, aes(x = fold, y = pcc, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Fold number", y = "PCC") +
  geom_hline(yintercept = 0.5, lty='dashed') +
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccBoxes_glm_scaled_perFold.pdf'),
       width = 6, height = 3.25)
g <- ggplot(res, aes(x = fold, y = rmse, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Fold number", y = "rmse") +
  theme_classic()  
ggsave(plot = g, file = paste0(OUT_DIR, 'rmseBoxes_glm_scaled_perFold.pdf'),
       width = 6, height = 3.25)

# Aggregate plots
# Per base position
tmp <- res[,.(fracAgree = length(which(pcc>=0.5))/.N),by = c('pos','modelType','trainType')]
g <- ggplot(tmp, aes(x = pos, y = fracAgree, fill = modelType)) +
  geom_bar(color = 'black', position = 'dodge', stat = 'identity') + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "Fraction of coumns in agreement") +
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccFracAgree_glm_scaled_perPos.pdf'),
       width = 6, height = 3.25)
tmp <- res[,.(fracAgree = length(which(pcc>=0.5))/.N),by = c('fold','modelType','trainType')]
g <- ggplot(tmp, aes(x = fold, y = fracAgree, fill = modelType)) +
  geom_bar(color = 'black', position = 'dodge', stat = 'identity') + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Fold number", y = "Fraction of coumns in agreement") +
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccFracAgree_glm_scaled_perFold.pdf'),
       width = 6, height = 3.25)

# Per fold
g <- ggplot(res, aes(x = fold, y = pcc, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "PCC") +
  geom_hline(yintercept = 0.5, lty='dashed') +
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccBoxes_glm_scaled_perFold.pdf'),
       width = 6, height = 3.25)
g <- ggplot(res, aes(x = fold, y = rmse, fill = modelType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "rmse") +
  theme_classic()  
ggsave(plot = g, file = paste0(OUT_DIR, 'rmseBoxes_glm_scaled_perFold.pdf'),
       width = 6, height = 3.25)

# Compare the RFs (scaled only) trained on the scaled vs. unscaled GLM alignments
res$glmType <- 'glm_scaled'
res.unscaled$glmType <- 'glm_unscaled'
res.cmp <- rbind(res, res.unscaled)
res.cmp$glmType <- factor(res.cmp$glmType)

# Fraction of examples correct
tmp <- res.cmp[trainType == 'rf_scaled',
               .(fracAgree = length(which(pcc>=0.5))/.N),
               by = c('pos','modelType','glmType')]
g <- ggplot(tmp, aes(x = pos, y = fracAgree, fill = glmType)) +
  geom_bar(color = 'black', position = 'dodge', stat = 'identity') + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "Fraction of coumns in agreement") +
  facet_wrap(~modelType)+
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccFracAgree_glm_scaledVsUnscaled_perPos.pdf'),
       width = 8, height = 3.25)
tmp <- res.cmp[trainType == 'rf_scaled',
               .(fracAgree = length(which(pcc>=0.5))/.N),
               by = c('fold','modelType','glmType')]
g <- ggplot(tmp, aes(x = fold, y = fracAgree, fill = glmType)) +
  geom_bar(color = 'black', position = 'dodge', stat = 'identity') + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  facet_wrap(~modelType)+
  labs(x = "Fold number", y = "Fraction of coumns in agreement") +
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccFracAgree_glm_scaledVsUnscaled_perFold.pdf'),
       width = 8, height = 3.25)

# PCC boxplots
g <- ggplot(res.cmp, aes(x = pos, y = pcc, fill = glmType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "PCC") +
  facet_wrap(~modelType) +
  geom_hline(yintercept = 0.5, lty='dashed') +
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccBoxes_glm_scaledVsUnscaled_perPos.pdf'),
       width = 8, height = 3.25)
g <- ggplot(res.cmp, aes(x = fold, y = pcc, fill = glmType)) +
  geom_boxplot(position = 'dodge', outlier.size = 0.1) + 
  scale_fill_brewer("Model type", palette = 'Dark2') + 
  labs(x = "Base position", y = "PCC") +
  facet_wrap(~modelType) +
  geom_hline(yintercept = 0.5, lty='dashed') +
  theme_classic()
ggsave(plot = g, file = paste0(OUT_DIR, 'pccBoxes_glm_scaledVsUnscaled_perFold.pdf'),
       width = 8, height = 3.25)
