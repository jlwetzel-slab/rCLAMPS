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
specs.protInfo <- specs[,c('prot','motif','zfNum','core','helix'),with=FALSE]
specs.protInfo <- specs.protInfo[!duplicated(specs.protInfo)]
specs.protInfo$motifNameStem.ffs <- sapply(specs.protInfo$prot, function(x) paste(x,'SOLEXA',sep = '_'))
specs.protInfo <- specs.protInfo[order(prot, zfNum)]
write.table(specs.protInfo, file = 'enuameh_perFinger_processedProtInfo.txt', 
            sep = '\t',quote = FALSE, row.names = FALSE)

# Output the per-finger PWMs as a table with the maximum finger first
# and removing the overlapping columns for adjacent fingers to align
# directly to the full PWMs (i.e. from flyFactor_dataset_A.txt)
specs.revFingerOrder <- specs[order(prot,-zfNum,motifPos)]
specs.revFingerOrder <- specs.revFingerOrder[,c("prot","motif","motifPos","strand","A","C","G","T"),with=FALSE]
specs.revFingerOrder <- specs.revFingerOrder[!duplicated(specs.revFingerOrder)]
write.table(specs.revFingerOrder, file = 'enuameh_perFinger_PWMs_reverseFingerOrder.txt', 
            sep = '\t',quote = FALSE, row.names = FALSE)

# Record the correct start/orientation information for the FFS proteins
#### NOTE:  There are inconsistencies in the strand given by Supplemental Dataset 1
####        vs. the correct strand in the corresponding PWM from flyfactor_dataset_A.txt, so 
####        only the start position is useful here currently.  Thus we directly align
####        each pair of PWMs to get the correct orientations.
startInfo <- specs[,.(start = min(motifPos) - 1, rev = ifelse(strand == -1, 1, 0)), by = 'prot']
startInfo <- startInfo[!duplicated(startInfo)]
write.table(startInfo, file = 'enuameh_perFinger_processedProt_startPosInfo.txt', 
            sep = '\t',quote = FALSE, row.names = FALSE)
