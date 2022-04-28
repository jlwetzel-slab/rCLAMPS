# Gather clustered ZFs based on the output of findDomains.py and subset 
# tables appropriately

library(data.table)

rm(list = ls())

inDoms <- fread('prot_seq_fewZFs_hmmerOut.txt')
inMotifs <- fread('motifTable_mostRecent_fewZFs.txt')

# Remove clusters of length less than 2
clustInfo <- inDoms[,.(clustLen = .N),by = c('prot','clustNum')]
doms <- merge(inDoms, clustInfo[clustLen >= 2 & clustLen <= 4], by = c('prot','clustNum'))

# Remove proteins with more than one cluster
nClusts <- doms[,.(nClusts = length(unique(clustNum))),by = c('prot')]
doms <- merge(doms, nClusts[nClusts == 1], by = c('prot'))

# Remove proteins for which we have no motif
doms <- doms[prot %in% inMotifs$TF_ID]
motifs <- inMotifs[TF_ID %in% doms$prot]

# Append the appropriate motif ID and TF name to the domain table
doms <- merge(doms, motifs[,c('TF_ID','Motif_ID','TF_Name','TF_Species'),with=FALSE], 
              by.x = 'prot', by.y = 'TF_ID', all.x = TRUE)
doms <- doms[order(prot, targStart)]

# Keep the subsetted domain and motif tables for organizing inputs to gibbs_GLM.py
write.table(doms, 'prot_seq_fewZFs_hmmerOut_clusteredOnly.txt',sep = '\t', quote = FALSE,row.names = FALSE)
write.table(motifs, 'motifTable_mostRecent_fewZFs_clusteredOnly.txt',sep = '\t', quote = FALSE,row.names = FALSE)

# Sanity check
print(length(unique(doms$prot)) == length(intersect(unique(motifs$TF_ID), unique(doms$prot))))
print(length(unique(motifs$TF_ID)) == length(intersect(unique(motifs$TF_ID), unique(doms$prot))))

# Overall properties of dataset
print(paste("There are",length(unique(doms$prot)),"distinct C2H2-ZFs proteins with a cis-bp",
            "motif and a single cluster of between 2 to 4 ZFs linked by at most 10 AAs b/w ZFs."))
print(paste("Thes proteins span:"))
print(paste(length(table(motifs$TF_Species)), "distinct species across", length(table(motifs$MSource_Identifier)), "distinct sources."))
print(paste("The distribution of cluster lengths is:"))
print(table(doms[,.(clustLen=.N), by = 'prot']$clustLen))

