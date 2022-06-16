# Create a contact map base on the zf-C2H2 structural alignment files

library(data.table)
rm(list = ls())

# Read in the aligned contact files and subset to base-contacts only (no backbone)
# Consider positions that contact any base in >= 10% of weighted structures
allBpos <- fread('zf-C2H2_weightedSum_distCut3.6_unionBases.txt')
allBpos <- allBpos[bType == 'base' & wFrac >= .1]  

# Read in the aligned contact files and subset to base-contacts only (no backbone)
# Consider all pairs that contact in >= 10% of weighted structures
pairs <- fread('zf-C2H2_weightedSum_distCut3.6.txt')
pairs <- pairs[cType == 'base' & frac >= 0.1 & mState %in% allBpos$mState] 
pairs <- pairs[bPos > 0]
pairs$bPos <- pairs$bPos - 1

#write.table(pairs)
