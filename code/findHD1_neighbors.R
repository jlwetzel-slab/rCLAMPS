# Explore the amino acid composition of sequences from across the various
# datasets and how many examples we could have for predicting specificities
# for proteins whose core sequence is HD-1 from another protein in the training set

library(data.table)
library(ggplot2)
library(RColorBrewer)
rm(list = ls())

# The core sequence positions being considered
EXCLUDE_CORE_OLAPS <- TRUE
if (EXCLUDE_CORE_OLAPS) {
  PLOT_DIR <- '../hd1-explore/0_excludeOlaps/'
} else {
  PLOT_DIR <- '../hd1-explore/'
}
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
  for (i in 1:nrow(matchTabs[dset == dset1])) {
    p <- matchTabs[dset == dset1]$prot[i]
    mutCore <- as.vector(matchTabs[dset == dset1][i,CORE_POS,with=FALSE])
    protOrig <- c()
    protMut <- c()
    to <- c()
    fr <- c()
    pos <- c()
    found <- FALSE
    for (j in 1:nrow(matchTabs[dset == dset2])) {
      origCore <- as.vector(matchTabs[dset == dset2][j,CORE_POS,with=FALSE])
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
        protMut <- c(protMut, matchTabs[dset == dset2]$prot[j])
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
  # and how many observations we have for each AA in the mutation
  # pair for that position
  
  d.perMut <- d[,.(nMuts = .N),by = c('corepos','from','to')]
  d.perMut <- merge(d.perMut, m, 
                    by.x = c('corepos','from'), by.y = c('position','aa'))
  d.perMut[,'dset':=NULL][]
  names(d.perMut)[length(d.perMut)] <- 'nObs.from'
  d.perMut <- merge(d.perMut, m, 
                    by.x = c('corepos','to'), by.y = c('position','aa'))
  d.perMut[,'dset':=NULL][]
  names(d.perMut)[length(d.perMut)] <- 'nObs.to'
  d.perMut$nObs.ge5 <- d.perMut$nObs.from >= 5 & d.perMut$nObs.to >= 5
  d.perMut$nObs.ge10 <- d.perMut$nObs.from >= 10 & d.perMut$nObs.to >= 10
  d.perMut$from <- factor(d.perMut$from, levels = AMINO)
  d.perMut$to <- factor(d.perMut$to, levels = AMINO)
  d.perMut$corepos <- factor(d.perMut$corepos, levels = CORE_POS)
  d.perMut
}

getFullHD1protsBarrera <- function(matchTabs) {
  coreMuts.barrera.exactProt <- NULL
  for (i in 1:nrow(matchTabs[dset == 'Barrera2016'])) {
    p <- matchTabs[dset == 'Barrera2016']$prot[i]
    pSub <- strsplit(p, split = '_')[[1]][1]
    subPos <- strsplit(p, split = '_')[[1]][2]
    origCore <- 
      as.vector(matchTabs[dset == 'CIS-BP' & prot == pSub,CORE_POS,with=FALSE])
    #print(paste(p, pSub))
    mutCore <- as.vector(matchTabs[dset == 'Barrera2016'][i,CORE_POS,with=FALSE])
    diff <- ifelse(origCore == mutCore, 0, 1)
    hDist <- sum(ifelse(origCore == mutCore, 0, 1))
    #print(origCore)
    if (hDist == 1) {
      thisPos <- which(diff == 1)
      fr <- as.character(origCore)[thisPos]; to <- as.character(mutCore)[thisPos]
      pos <- CORE_POS[thisPos]
      coreMuts.barrera.exactProt <- rbind(coreMuts.barrera.exactProt, 
                                          data.table(prot = pSub, mutProt = p,
                                                     corepos = pos, from = fr, to = to))
    }
  }
  coreMuts.barrera.exactProt
}

matchTabs <- NULL
for (i in 1:length(MATCHFILES)) {
  tmp <- fread(MATCHFILES[i])
  tmp$dset <- GRPS[i]
  names(tmp) <- c('prot', paste0('A.',2:(ncol(tmp)-1)),'dset')
  matchTabs <- rbind(matchTabs, tmp)
}
matchTabs <- matchTabs[,c('prot',CORE_POS,'dset'),with = FALSE]

if (EXCLUDE_CORE_OLAPS) {
  matchTabs$coreSeq <- sapply(1:nrow(matchTabs), function(i) {
    paste0(as.character(matchTabs[i,CORE_POS,with=FALSE]),collapse = '')
  })
  testCores <- unique(matchTabs[dset == 'Barrera2016' | dset == 'Chu2012']$coreSeq)
  trainCores <- unique(matchTabs[dset == 'CIS-BP' & !(coreSeq %in% testCores)]$coreSeq)
  trainSetSz.orig <- nrow(matchTabs[dset == 'CIS-BP'])
  matchTabs <- rbind(matchTabs[dset == 'CIS-BP' & coreSeq %in% trainCores],
                     matchTabs[dset %in% c('Barrera2016','Chu2012')])
  matchTabs[,coreSeq:=NULL][]
}
trainSetSz.final <- nrow(matchTabs[dset == 'CIS-BP'])
print(paste("Removed",trainSetSz.orig-trainSetSz.final,
            "overlapping proteins from the training set."))

# How many times do we observe each AA in each position the two 
# primary datasets?
matchTabs.melted <- NULL
for (ds in c('CIS-BP','Chu2012')) {
  tmp <- matchTabs[dset == ds,c('prot',CORE_POS),with=FALSE]
  tmp <- melt(tmp, id.vars = 'prot', measure.vars = CORE_POS,
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
  matchTabs.melted[,.(nObs = .N),by = c('dset','position','aa')]

# How much does adding the Chu2012 dataset get us w.r.t 47, 50 and 54?
g <- ggplot(matchTabs.melted[position %in% paste0('A.',c('47','50','54'))], aes(x = aa)) + 
  geom_bar(color = 'black', aes(fill = dset)) + 
  facet_wrap(~position, scales = 'free_y', nrow = 1) + 
  scale_fill_brewer("",palette = 'Dark2') +
  labs(x = 'Amino acid') +
  theme_classic() + 
  theme(legend.position = 'top')
ggsave(plot = g, file = paste0(PLOT_DIR, 'aa_perPos_both.pdf'),
       height = 3.6, width = 9)

# How many proteins in the Barerra dataset vary from their corresponding native protein
# in one of the core sequence positions?
fname <- paste0(PLOT_DIR,'0_hd1fullProts_Barrera2016_CIS-BP.txt')
if (!file.exists(fname)) {
  coreMuts.barrera.exactProt <- getFullHD1protsBarrera(matchTabs)
  write.table(coreMuts.barrera.exactProt, file = fname,
              quote = FALSE, row.names = FALSE, sep = '\t')
} else {
  coreMuts.barrera.exactProt <- fread(fname)
}

# How many core sequences in the Barrera dataset are HD1 from anything 
# in the CIS-BP dataset?
fname <- paste0(PLOT_DIR,'0_hd1cores_Barrera2016_CIS-BP.txt')
if (!file.exists(fname)) {
  coreMuts.barrera <- searchForHD1(matchTabs,'Barrera2016','CIS-BP')
  write.table(coreMuts.barrera, file = fname,
              quote = FALSE, row.names = FALSE, sep = '\t')
} else {
  coreMuts.barrera <- fread(fname)
}

# How many core sequences in the Chu2012 dataset are HD1 from something 
# in the CIS-BP dataset?
fname <- paste0(PLOT_DIR,'0_hd1cores_Chu2012_CIS-BP.txt')
if (!file.exists(fname)) {
  coreMuts.chu2012 <- searchForHD1(matchTabs,'Chu2012','CIS-BP')
  write.table(coreMuts.chu2012, file = fname,
              quote = FALSE, row.names = FALSE, sep = '\t')
} else {
  coreMuts.chu2012 <- fread(fname)
}

# How many core sequences in the Chu2012 dataset are HD1 from something 
# else also in the Chu2012 dataset?
fname <- paste0(PLOT_DIR,'0_hd1cores_Chu2012_Chu2012.txt')
if (!file.exists(fname)) {
  coreMuts.chu2012.self <- searchForHD1(matchTabs,'Chu2012','Chu2012')
  write.table(coreMuts.chu2012.self, file = fname,
              quote = FALSE, row.names = FALSE, sep = '\t')
} else {
  coreMuts.chu2012.self <- fread(fname)
}

# Summarize the number of times each HD-1 mutation is observed between 
# any core seq in either Barrera (or Chu) and any core seq in CIS-BP 
# and ask how much power whether we have at least 5 (or 10) mutations
# for both of the corresponding amino acids in the corresponding core position.
mutSummary.barrera.cbp <- 
  getMutSummaryTab(coreMuts.barrera,matchTabs.melted.summary[dset == 'CIS-BP'])
mutSummary.chu2012.cbp <- 
  getMutSummaryTab(coreMuts.chu2012,matchTabs.melted.summary[dset == 'CIS-BP'])

g <- ggplot(mutSummary.barrera.cbp, aes(x = to, y = from)) + 
  geom_tile(color = 'black', aes(fill = nMuts)) + 
  geom_point(data = mutSummary.barrera.cbp[nObs.ge5 == TRUE], size = 0.5) +
  labs(y = "Amino acid (from)", x ="Amino acid (to)") +
  scale_fill_gradient2("Number of\nmutation\nexamples", low = 'blue', high = 'red')+
  facet_wrap(~corepos) +
  theme_classic()
ggsave(plot = g, width = 8, height = 8,
       file = paste0(PLOT_DIR, 'nMutsHD-1_Barrera2016_CIS-BP_nObsge5.pdf'))
g <- ggplot(mutSummary.barrera.cbp, aes(x = to, y = from)) + 
  geom_tile(color = 'black', aes(fill = nMuts)) + 
  geom_point(data = mutSummary.barrera.cbp[nObs.ge10 == TRUE], size = 0.5) +
  labs(y = "Amino acid (from)", x ="Amino acid (to)") +
  scale_fill_gradient2("Number of\nmutation\nexamples", low = 'blue', high = 'red')+
  facet_wrap(~corepos) +
  theme_classic()
ggsave(plot = g, width = 8, height = 8,
       file = paste0(PLOT_DIR, 'nMutsHD-1_Barrera2016_CIS-BP_nObsge10.pdf'))
g <- ggplot(mutSummary.chu2012.cbp, aes(x = to, y = from)) + 
  geom_tile(color = 'black', aes(fill = nMuts)) + 
  geom_point(data = mutSummary.chu2012.cbp[nObs.ge5 == TRUE], size = 0.5) +
  labs(y = "Amino acid (from)", x ="Amino acid (to)") +
  scale_fill_gradient2("Number of\nmutation\nexamples", low = 'blue', high = 'red')+
  facet_wrap(~corepos) +
  theme_classic()
ggsave(plot = g, width = 5, height = 4,
       file = paste0(PLOT_DIR, 'nMutsHD-1_Chu2012_CIS-BP_nObsge5.pdf'))
g <- ggplot(mutSummary.chu2012.cbp, aes(x = to, y = from)) + 
  geom_tile(color = 'black', aes(fill = nMuts)) + 
  geom_point(data = mutSummary.chu2012.cbp[nObs.ge10 == TRUE], size = 0.5) +
  labs(y = "Amino acid (from)", x ="Amino acid (to)") +
  scale_fill_gradient2("Number of\nmutation\nexamples", low = 'blue', high = 'red')+
  facet_wrap(~corepos) +
  theme_classic()
ggsave(plot = g, width = 5, height = 4,
       file = paste0(PLOT_DIR, 'nMutsHD-1_Chu2012_CIS-BP_nObsge10.pdf'))

# What about predicting HD1-core mutations found in Chu using other HD1 proteins
# found in Chu et al., but with model training on CIS-BP
mutSummary.chu2012self.cbp <- 
  getMutSummaryTab(coreMuts.chu2012.self,matchTabs.melted.summary[dset == 'CIS-BP'])

g <- ggplot(mutSummary.chu2012self.cbp, aes(x = to, y = from)) + 
  geom_tile(color = 'black', aes(fill = nMuts)) + 
  geom_point(data = mutSummary.chu2012self.cbp[nObs.ge5 == TRUE], size = 0.5) +
  labs(y = "Amino acid (from)", x ="Amino acid (to)") +
  scale_fill_gradient2("Number of\nmutation\nexamples", low = 'blue', high = 'red')+
  facet_wrap(~corepos) +
  theme_classic()
ggsave(plot = g, width = 8, height = 3,
       file = paste0(PLOT_DIR, 'nMutsHD-1_Chu2012self_CIS-BP_nObsge5.pdf'))
g <- ggplot(mutSummary.chu2012self.cbp, aes(x = to, y = from)) + 
  geom_tile(color = 'black', aes(fill = nMuts)) + 
  geom_point(data = mutSummary.chu2012self.cbp[nObs.ge10 == TRUE], size = 0.5) +
  labs(y = "Amino acid (from)", x ="Amino acid (to)") +
  scale_fill_gradient2("Number of\nmutation\nexamples", low = 'blue', high = 'red')+
  facet_wrap(~corepos) +
  theme_classic()
ggsave(plot = g, width = 8, height = 3,
       file = paste0(PLOT_DIR, 'nMutsHD-1_Chu2012self_CIS-BP_nObsge10.pdf'))

