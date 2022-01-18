# Explore the amino acid composition of sequences from across the various
# datasets and how many examples we could have for predicting specificities
# for proteins whose core sequence is HD-1 from another protein in the training set
# if we were to randomly split the Chu dataset and include 1/2 of it in our 
# training set.

library(data.table)
library(ggplot2)
library(RColorBrewer)
rm(list = ls())

RAND_SEED <- 23646326
PLOT_DIR <- '../hd1-explore/0_splitChu2012-trainTest/'
# The core sequence positions being considered
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
CORE_POS <- paste0('A.',c('2','3','4','5','47','50','51','54','55'))
GRPS <- c('CIS-BP', 'Barrera2016', 'Chu2012')
MATCHFILES <- c('../cis_bp/homeodomains_cisbp_hasPWM.matchTab.txt',
                '../cis_bp/homeodomains_barreraMuts_hasPWM.matchTab.txt',
                '../Chu2012_bsSelections/allProts_hmmer_out__hasPWM.matchTab.txt')
AMINO <- c('A','C','D','E','F','G','H','I','K','L',
           'M','N','P','Q','R','S','T','V','W','Y')

searchForHD1 <- function(matchTabs, dset1, dset2) {
  hd1Info <- NULL
  for (i in 1:nrow(matchTabs[grp == dset1])) {
    p <- matchTabs[grp == dset1]$prot[i]
    mutCore <- as.vector(matchTabs[grp == dset1][i,CORE_POS,with=FALSE])
    protOrig <- c()
    protMut <- c()
    to <- c()
    fr <- c()
    pos <- c()
    found <- FALSE
    for (j in 1:nrow(matchTabs[grp == dset2])) {
      origCore <- as.vector(matchTabs[grp == dset2][j,CORE_POS,with=FALSE])
      diff <- ifelse(origCore == mutCore, 0, 1)
      #print(diff)
      hDist <- sum(ifelse(origCore == mutCore, 0, 1))
      if (hDist == 1) {
        #print("Found one!")
        found <- TRUE
        thisPos <- which(diff == 1)
        this.fr <- as.character(origCore)[thisPos]; 
        this.to <- as.character(mutCore)[thisPos]
        this.pos <- CORE_POS[thisPos]
        protOrig <- c(protOrig, p)
        protMut <- c(protMut, matchTabs[grp == dset2]$prot[j])
        to <- c(to,this.to); 
        fr <- c(fr,this.fr); 
        pos <- c(pos,this.pos)
      }
    }
    if (found) {
      tmp <- data.table(prot = protOrig, mutProt = protMut, corepos = pos, 
                        from = fr, to = to)
      hd1Info <- rbind(hd1Info, tmp)
    }
    print(paste(i, found, nrow(hd1Info)))
  }
  hd1Info
}

getMutSummaryTab <- function(d, m) {
  # Takes as input an HD-1 mutation data table (d) in melted format
  # and a match table summary (m) and returns a table summarizing
  # the number of observations available for each mutation type 
  # and how many observations we have for the AA being substituted to
  
  d.perMut <- d[,.(nMuts = .N),by = c('corepos','from','to')]
  d.perMut <- merge(d.perMut, m, 
                    by.x = c('corepos','to'), by.y = c('position','aa'))
  d.perMut[,c('dset'):=NULL][]
  names(d.perMut)[length(d.perMut)] <- 'nObs.to'
  d.perMut$nObs.ge5 <- d.perMut$nObs.to >= 5
  d.perMut$nObs.ge10 <- d.perMut$nObs.to >= 10
  d.perMut$from <- factor(d.perMut$from, levels = AMINO)
  d.perMut$to <- factor(d.perMut$to, levels = AMINO)
  d.perMut$corepos <- factor(d.perMut$corepos, levels = CORE_POS)
  d.perMut
}

# Seed the RNG
set.seed(RAND_SEED)

matchTabs <- NULL
for (i in 1:length(MATCHFILES)) {
  tmp <- fread(MATCHFILES[i])
  tmp$dset <- GRPS[i]
  names(tmp) <- c('prot', paste0('A.',2:(ncol(tmp)-1)),'dset')
  matchTabs <- rbind(matchTabs, tmp)
}
matchTabs <- matchTabs[,c('prot',CORE_POS,'dset'),with = FALSE]

# Break up into training and test sets.  We start with all CIS-BP proteins
# and a random subset of proteins corresponding to 1/2 of the distinct core 
# sequences that occur in the Chu 2012 dataset as our training set.  We then 
# remove any proteins whose core sequences overlap with the Barrera 2016 proteins
# from this training set.  Our test set will then be the entire Barrera 2016 
# dataset plus the proteins corresponding to the remaining 1/2 of distinct core
# sequences from the Chu 2012 dataset.
matchTabs$coreSeq <- sapply(1:nrow(matchTabs), function(i) {
  paste0(as.character(matchTabs[i,CORE_POS,with=FALSE]),collapse = '')
})
chuCores <- unique(matchTabs[dset == 'Chu2012']$coreSeq)
chu.train <- sample(chuCores, size = length(unique(chuCores))/2)
chu.test <- setdiff(chuCores, chu.train)
testCores <- union(chu.test, unique(matchTabs[dset == 'Barrera2016']$coreSeq))
#trainCores <- setdiff(unique(matchTabs[dset == 'CIS-BP']$coreSeq), testCores)
matchTabs$grp <- sapply(matchTabs$coreSeq, function(c) {
  if (c %in% testCores) 'test' else 'train'
})
print("Distribution of training and testing examples across datasetst:")
print(table(matchTabs$dset, matchTabs$grp))
write.table(matchTabs[grp == 'test']$prot, file = paste0(PLOT_DIR, 'testProts.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(matchTabs, file = paste0(PLOT_DIR, '0_matchTabs.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)

# How many times do we observe each AA in each position the training set?
matchTabs.melted <- NULL
for (ds in c('CIS-BP','Chu2012')) {
  tmp <- matchTabs[dset == ds & grp == 'train',c('prot','grp',CORE_POS),with=FALSE]
  tmp <- melt(tmp, id.vars = c('prot','grp'), measure.vars = CORE_POS,
              variable.name = 'position', value.name = 'aa')
  tmp$dset <- ds
  matchTabs.melted <- rbind(matchTabs.melted,tmp)
  rm(tmp)
  
  g <- ggplot(matchTabs.melted[dset == ds], aes(x = aa)) + 
    geom_bar(color = 'black') + 
    facet_wrap(~position, scales = 'free_y') + 
    labs(x = 'Amino acid') +
    theme_classic()
  ggsave(plot = g, file = paste0(PLOT_DIR, 'aa_perPos_',ds,'.pdf'),
         height = 6, width = 8)
}
matchTabs.melted.summary <- 
  matchTabs.melted[,.(nObs = .N),by = c('dset','grp','position','aa')]

# How much does adding the Chu2012 dataset get us w.r.t 47, 50 and 54?
g <- ggplot(matchTabs.melted[position %in% paste0('A.',c('47','50','54'))], 
            aes(x = aa)) + 
  geom_bar(color = 'black', aes(fill = dset)) + 
  facet_wrap(~position, scales = 'free_y', nrow = 1) + 
  scale_fill_brewer("",palette = 'Dark2') +
  labs(x = 'Amino acid') +
  theme_classic() + 
  theme(legend.position = 'top')
ggsave(plot = g, file = paste0(PLOT_DIR, 'aa_perPos_both.pdf'),
       height = 3.6, width = 9)

# How many core sequences in the test dataset are HD1 from anything 
# in the training dataset?
fname <- paste0(PLOT_DIR,'0_hd1cores_test_train.txt')
if (!file.exists(fname)) {
  coreMuts.test <- searchForHD1(matchTabs,'test','train')
  write.table(coreMuts.test, file = fname,
              quote = FALSE, row.names = FALSE, sep = '\t')
} else {
  coreMuts.test <- fread(fname)
}

# Summarize the number of times each HD-1 mutation is observed between 
# any core seq in the training set and a core sequence in the test, 
# and ask whether we have at least 5 (or 10) examples of the amino acid 
# being substituted to in the traing set

# for both of the corresponding amino acids in the corresponding core position.
mutSummary <- 
  getMutSummaryTab(coreMuts.test,matchTabs.melted.summary[grp == 'train'])

g <- ggplot(mutSummary, aes(x = to, y = from)) + 
  geom_tile(color = 'black', aes(fill = nMuts)) + 
  geom_point(data = mutSummary[nObs.ge5 == TRUE], size = 0.5) +
  labs(y = "Amino acid (from)", x ="Amino acid (to)") +
  scale_fill_gradient2("Number of\nmutation\nexamples", low = 'blue', high = 'red')+
  facet_wrap(~corepos) +
  theme_classic()
ggsave(plot = g, width = 8, height = 8,
       file = paste0(PLOT_DIR, 'nMutsHD-1_test_train_nObsge5.pdf'))
g <- ggplot(mutSummary, aes(x = to, y = from)) + 
  geom_tile(color = 'black', aes(fill = nMuts)) + 
  geom_point(data = mutSummary[nObs.ge10 == TRUE], size = 0.5) +
  labs(y = "Amino acid (from)", x ="Amino acid (to)") +
  scale_fill_gradient2("Number of\nmutation\nexamples", low = 'blue', high = 'red')+
  facet_wrap(~corepos) +
  theme_classic()
ggsave(plot = g, width = 8, height = 8,
       file = paste0(PLOT_DIR, 'nMutsHD-1_test_train_nObsge10.pdf'))