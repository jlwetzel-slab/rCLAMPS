# Examine the fly factor survey data

library(data.table)

specs <- fread('enuameh_perFinger.txt') 
specs <- specs[order(prot,motifPos)]
table(specs$zfNum)

# Can see based on comparison of these motif positions and 
# the motifs in 'flyfactor_dataset_A.txt' (from FlyFactorSurvey website)
# the these flyfactor PWMs are 1-indexed (while motifs are 0 indexed in my code) and 
# that these are roughly +1 smoothed versions of the corresponing columns in 'flyfactor_dataset_A.txt'
# e.g.,
motifPos  <- as.numeric(specs[prot == 'Blimp-1' & motifPos == 3,c('A','C','G','T'),with = FALSE])
# Compare to 7th position of >CG10267_SOLEXA_5_FBgn0037446 in 'flyfactor_dataset_A.txt'
all(abs((c(551,562,245,581)+1)/sum((c(551,562,245,581)+1)) - motifPos) < 1e-3)

# Remove proteins if their annotated ZFs start before position 1 of the 'flyfactor_dataset_A.txt' motif
keepProts <- specs[,.(keepProt = all(motifPos>0)), by = 'prot']
keepProts <- keepProts[keepProt == TRUE]$prot
specs <- specs[prot %in% keepProts]

# Order based on finger positions and output table to read into gibbs sampler
specs.protInfo$motifNameStem.ffs <- sapply(specs.protInfo$prot, function(x) paste(x,'SOLEXA',sep = '_'))
write.table(specs.protInfo, file = 'enuameh_perFinger_processedProtInfo.txt', 
            sep = '\t',quote = FALSE, row.names = FALSE)

# Record the correct start/orientation information for the FFS proteins
startInfo <- specs[,.(start = min(motifPos) - 1, rev = ifelse(strand == -1, 1, 0)), by = 'prot']
startInfo <- startInfo[!duplicated(startInfo)]
write.table(specs.protInfo, file = 'enuameh_perFinger_processedProt_startPosInfo.txt', 
            sep = '\t',quote = FALSE, row.names = FALSE)
