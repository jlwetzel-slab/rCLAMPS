# Gather clustered ZFs based on the output of findDomains.py and subset 
# tables appropriately

library(data.table)
library(ggplot2)

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

# How long are the motifs relative to the cluster length?
doms <- merge(doms, fread('PWM_lengths.txt'), by = 'Motif_ID')
doms <- doms[order(prot, targStart)]

# Sanity check
print(length(unique(doms$prot)) == length(intersect(unique(motifs$TF_ID), unique(doms$prot))))
print(length(unique(motifs$TF_ID)) == length(intersect(unique(motifs$TF_ID), unique(doms$prot))))

# Keep the subsetted domain and motif tables for organizing inputs to gibbs_GLM.py
write.table(doms, 'prot_seq_fewZFs_hmmerOut_clusteredOnly.txt',sep = '\t', quote = FALSE,row.names = FALSE)
write.table(motifs, 'motifTable_mostRecent_fewZFs_clusteredOnly.txt',sep = '\t', quote = FALSE,row.names = FALSE)

# Overall properties of dataset
print(paste("There are",length(unique(doms$prot)),"distinct C2H2-ZFs proteins with a cis-bp",
            "motif and a single cluster of between 2 to 4 ZFs linked by at most 10 AAs b/w ZFs."))
print(paste("These proteins span:"))
print(paste(length(unique(doms$coreSeq)), "distinct base-contacting amino acid combinations."))
print(paste(length(table(motifs$TF_Species)), "distinct species across", length(table(motifs$MSource_Identifier)), "distinct sources."))
print(paste("The distribution of cluster lengths is:"))
print(table(doms[,.(clustLen=.N), by = 'prot']$clustLen))

# Plot the distribution of (motifLen-1)/clustLen
g <- ggplot(doms[,.SD[1], by = 'prot'], aes(x = clustLen, y = motifLen, group = clustLen)) + 
  geom_boxplot(outlier.colour = 'white') + 
  geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4) + 
  labs(x = "Number of ZFs in cluster", "Motif length") + 
  theme_classic()
ggsave(plot = g, file = '0_nZFs_V_motifLen.pdf', width = 5, height = 5)

# Remove motif if not greater than 3*clustLen in length
doms <- doms[motifLen > 3*clustLen]
motifs <- motifs[TF_ID %in% doms$prot]

# Sanity check
print(length(unique(doms$prot)) == length(intersect(unique(motifs$TF_ID), unique(doms$prot))))
print(length(unique(motifs$TF_ID)) == length(intersect(unique(motifs$TF_ID), unique(doms$prot))))

# Overall properties of dataset
print("Considering only motifs where motifLen > 3*clustLen:")
print(paste("There are",length(unique(doms$prot)),"distinct C2H2-ZFs proteins with a cis-bp",
            "motif and a single cluster of between 2 to 4 ZFs linked by at most 10 AAs b/w ZFs."))
print(paste("These proteins span:"))
print(paste(length(unique(doms$coreSeq)), "distinct base-contacting amino acid combinations."))
print(paste(length(table(motifs$TF_Species)), "distinct species across", length(table(motifs$MSource_Identifier)), "distinct sources."))
print(paste("The distribution of cluster lengths is:"))
print(table(doms[,.(clustLen=.N), by = 'prot']$clustLen))

# Plot the distribution of (motifLen-1)/clustLen
g <- ggplot(doms[,.SD[1], by = 'prot'], aes(x = clustLen, y = motifLen, group = clustLen)) + 
  geom_boxplot(outlier.colour = 'white') + 
  geom_jitter(height = 0, width = 0.1, size = 1, alpha = 0.4) + 
  labs(x = "Number of ZFs in cluster", "Motif length") + 
  theme_classic()
ggsave(plot = g, file = '0_nZFs_V_motifLen_removeTooShort.pdf', width = 5, height = 5)

# Keep the subsetted domain and motif tables for organizing inputs to gibbs_GLM.py
arraySeqs <- doms[,.(arraySeq = paste(coreSeq, collapse = '')), by = 'prot']
doms <- merge(doms, arraySeqs, by = 'prot')
doms <- doms[order(prot, targStart)]
write.table(doms, 'prot_seq_fewZFs_hmmerOut_clusteredOnly_removeTooShort.txt',sep = '\t', quote = FALSE,row.names = FALSE)
write.table(motifs, 'motifTable_mostRecent_fewZFs_clusteredOnly_removeTooShort.txt',sep = '\t', quote = FALSE,row.names = FALSE)

# Which are the really long motifs??
longMotifs <- motifs[TF_ID %in% doms[motifLen > 4*clustLen]$prot]
longMotifs.doms <- doms[motifLen > 4*clustLen]

# How many distinct arrays are there (in terms of base-contacting positions)?
distinctArrays <- doms[,.(nProts = length(unique(prot))), by = 'arraySeq']
distinctArrays <- distinctArrays[order(-nProts)]
print(paste("There are", nrow(distinctArrays), "distinct ZF arrays (considering only base-contacting positions)."))

length(sort(table(names(sapply(doms$arraySeq, nchar)))))

# Which motifs have length exactly 3*nZFs + 1?


