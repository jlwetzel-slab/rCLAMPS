library(data.table)
library(ggplot2)
library(RColorBrewer)

rm(list = ls())

OUT_TAB <- '../cis_bp/'
PLOT_DIR <- '../cis_bp/plots/'

# From '../../../../../cisBP/gatherHomeoboxMotifs.R' 
motifInfo <- fread('../../../../../cisBP/hboxOnlyInfo_includeInVivo.txt')

# Remove JASPAR as they nearly all overlap with a primary source
# and remove Transfac since we don't have the lisence
motifInfo <- motifInfo[MSource_Identifier != 'Transfac']
motifInfo <- motifInfo[MSource_Identifier != 'JASPAR']

# Set aside motifs corresponding to mutational studies
mutants <- motifInfo[MSource_Identifier %in% c('Barrera2016')]
mutants <- mutants[!grep('_REF$', mutants$DBID.2)]

# Make a table where we keep only 1 motif per TF;
# Below are rules for which to keep from particular sources in order 
# to exclude perturbation studies; after removing motifs corresponding 
# to these, we simply keep the most recent motif for each distinct TF
sourceRules <- list('Barrera2016' = '_REF$', 'Yin2017'= 'eDBD_HT-SELEX',
                    'Yin2017b' = '_Cytosine$','FlyFactorSurvey' = '_Cell_')
for (sName in names(sourceRules)) {
  keep <- motifInfo[grep(sourceRules[[sName]], motifInfo$DBID.2)]
  motifInfo <- motifInfo[MSource_Identifier != sName]
  motifInfo <- rbind(motifInfo, keep)
}

# Prioritize in vitro data first (keeping most recent only)
inVivoTypes <- c('ChIP-seq','Dap-seq','Misc')
motifInfo.inVivo <- motifInfo[MSource_Type %in% inVivoTypes]
motifInfo <- motifInfo[!(MSource_Type %in% inVivoTypes)]
motifInfo <- motifInfo[order(motifInfo$MSource_Year, decreasing = TRUE)]
motifInfo <- unique(motifInfo, by = 'TF_Name')

# How many are covered in vivo but not in vitro?  Add back in 
# if not covered in vitro
setdiff(motifInfo.inVivo$TF_Name,motifInfo$TF_Name) #There are 27
motifInfo.inVivo <- 
  motifInfo.inVivo[TF_Name %in% setdiff(motifInfo.inVivo$TF_Name,motifInfo$TF_Name)]
motifInfo.inVivo <- motifInfo.inVivo[order(motifInfo.inVivo$MSource_Year, 
                                           decreasing = TRUE)]
motifInfo.inVivo <- unique(motifInfo.inVivo, by = 'TF_Name')
motifInfo <- rbind(motifInfo, motifInfo.inVivo)
write.table(motifInfo, file = paste0(OUT_TAB,'motifTable_mostRecent_noMuts.txt'),
            quote = FALSE, row.names = FALSE, sep = '\t')

# Make a table for the PBMs with mutants
write.table(mutants, file = paste0(OUT_TAB,'motifTable_Barrera2016_mutsOnly.txt'),
            quote = FALSE, row.names = FALSE, sep = '\t')
