# A python script to run MEME on all of the HD binding site selection 
# data in the supplemental file from Chu et al. 2012 Supplemental Table S6

import os
import numpy as np
import multiprocessing as mp

SUPP_TAB_FILE = '../Chu2012_SupTableS6.txt'
FASTA_DIR = '../fasta/'
MEME_OUT_DIR = '../MEME_outputs/'
MOTIF_SUMM_TAB = '../MEME_motifSummaryInfo.txt'
BG_SUMM_TAB = '../MEME_motif_bgInfo.txt'
MOTIF_TAB = '../MEME_allMotifs.txt'
MOTIF_TAB_RW = '../MEME_allMotifs_rw.txt'
MOTIF_TAB_RW_LOG2 = '../MEME_allMotifs_rw_log2.txt'

RUN_MEME = False #True#

def makeFastaFiles(supTabFile, fastaDir):
    # Creates a set of fasta files based on a tab-delimited text 
    # version of the supplemental table S6

    seqs = []
    fin = open(supTabFile, 'r')
    selectionLabs = fin.readline().strip().split('\t')
    nSelections = len(selectionLabs)
    fin.readline()
    fin.readline()
    fin.readline()
    print nSelections
    labels = [[] for i in range(nSelections)]
    seqs = [[] for i in range(nSelections)]
    for k, line in enumerate(fin):
        l = line.split('\t')
        l[len(l)-1] = l[len(l)-1].strip()
        #print k, len(l)
        if l[0] != '' and l[0][0] == '>':
            for i in range(nSelections):
                if l[i] != '':
                    labels[i].append(l[i])
        else:
            for i in range(nSelections):
                if l[i] != '':
                    seqs[i].append(l[i])  # Don't remove the NN flanks
                    #seqs[i].append(l[i][2:-2])  # Remove the NN flanks
    fin.close()
    maxSeqs = max([len(x) for x in seqs])
    
    for i in range(nSelections):
        fname = fastaDir+selectionLabs[i]+'.fasta'
        fout = open(fname, 'w')
        for j in range(len(labels[i])):
            fout.write('%s\n' %labels[i][j])
            fout.write('%s\n' %seqs[i][j])
        fout.close()

    return maxSeqs


def parseMemeTxt(fname):
    # Returns a dictionary of information extracted from
    # the meme.txt file at path fname

    def parseMotifInfo(fin, line):
        # Parse a motif-instance where fin has just read the
        # line labeled as the start of a motif entry

        l = line.strip().split()
        mInfo = {l[i]:l[i+2] for i in range(3,13,3)}
        mInfo['motifNum'] = l[1]

        # Get the set of alignment sequence info
        l = fin.readline().strip().split()
        while len(l) == 0 or l[0] != "Sequence":
            l = fin.readline().strip().split()
        l = [fin.readline().strip().split() for i in range(2)][1]
        mInfo['aliSeqs'] = []
        while l[0][0] != '-':
            mInfo['aliSeqs'].append('\t'.join(l[:3]+[l[5]]))
            l = fin.readline().strip().split()

        # Find the pfm info
        l = fin.readline().strip().split()
        while len(l) == 0 or l[0] != "letter-probability":
            l = fin.readline().strip().split()
        l = fin.readline().strip().split()
        mInfo['pfm'] = []
        while l[0][0] != '-':
            mInfo['pfm'].append([x for x in l])
            l = fin.readline().strip().split()
        return mInfo

    info = {}
    info['motifInfo'] = []  # List of motif info dictionaries
    fin = open(fname)

    # Get 0-order letter frequencies
    l = fin.readline().strip().split()
    while len(l) == 0 or l[0] != "Letter":
        l = fin.readline().strip().split()
    l = fin.readline().strip().split()
    info['letterFreq'] = {l[i]: l[i+1] for i in range(0,len(l),2)}
    l = [fin.readline().strip().split() for i in range(2)][1]
    info['bgFreq'] = {l[i]: l[i+1] for i in range(0,len(l),2)}

    # Get motifInfo for each motif in order of appearance
    line = fin.readline()
    while len(line) > 0:
        l = line.strip().split()
        if len(l) > 0 and l[0] == "MOTIF":
            info['motifInfo'].append(parseMotifInfo(fin, line))
        line = fin.readline()

    return info

def runMEME(fasta, outDir, minw = 4, maxw = 10, nmotifs = 1):
    # Run MEME with on DNA with zoops

    sysCall = 'meme %s -dna -mod zoops -nmotifs %d ' %(fasta, nmotifs) + \
        '-minw %d -maxw %d -revcomp -oc %s -brief 2000' %(minw, maxw, outDir)
    os.system(sysCall)

def makeWeightedMotif(aliInfo, norm = True, log2 = False):
    # Give the aliSeqsInfo returned by parseMemeFile, returns a 
    # reweighted motif in the same way as Chu et. al

    ### Note: 'N' is treated as 1/4 of each base

    BASES = ['A','C','G','T']
    mwid = len(aliInfo[0].split('\t')[3])
    m = np.zeros((mwid,4), dtype = 'float')
    ambiguous = 0
    for line in aliInfo:
        l = line.split('\t')
        wt, seq = int(l[0].split('.')[1]), l[3]
        for j in range(mwid):
            if seq[j] == 'N':
                if log2:
                    m[j] += np.log2(wt)
                else:
                    m[j] += wt
            else:
                if log2:
                    m[j,BASES.index(seq[j])] += np.log2(wt)
                else:
                    m[j,BASES.index(seq[j])] += wt
    if norm:
        for j in range(mwid):
            m[j] = m[j]/m[j].sum()


    return m

def main():

    # Create the fasta files out of the supplemental table
    if not os.path.exists(FASTA_DIR):
        os.makedirs(FASTA_DIR)
    if not os.path.exists(MEME_OUT_DIR):
        os.makedirs(MEME_OUT_DIR)
    maxSeqs = makeFastaFiles(SUPP_TAB_FILE, FASTA_DIR)


    # Run MEME in parallel on each of the fasta files
    fnames = sorted(os.listdir(FASTA_DIR))#[:1]
    #fname = fnames[0]
    #runMEME(FASTA_DIR+fname,MEME_OUT_DIR+fname.split('.')[0])
    
    #"""
    if RUN_MEME:
        pool = mp.Pool(processes = mp.cpu_count())
        results = \
            [pool.apply_async(runMEME, 
                              args = (FASTA_DIR+fname,
                                      MEME_OUT_DIR+fname.split('.')[0])) \
            for fname in fnames]
        # Force evaluation of the processes as they finish
        [r.get() for r in results]
    #"""


    # Parse the MEME output and place into tables, etc ...
    print "Parsing meme output ..."
    msout = open(MOTIF_SUMM_TAB, 'w')   # Motif summary file
    bsout = open(BG_SUMM_TAB, 'w')      # Nucleotide frequencies
    mout = open(MOTIF_TAB, 'w')         # Complete set of motifs
    mout_rw = open(MOTIF_TAB_RW, 'w')   # Complete set of motifs reweighted as ins Chu et al.
    mout_rw_log2 = open(MOTIF_TAB_RW_LOG2, 'w')  # Complete set of motifs reweighted with logs of counts
    mSummLabs = ['uniprot','motifNum','sites','width','llr','E-value']
    bSummLabs = ['uniprot'] + ['A','C','G','T'] + ['freqType']
    msout.write('\t'.join(['label','consensus','sites','width','llr','E-value']) + '\n')
    bsout.write('\t'.join(['label','A','C','G','T','freqType']) + '\n')
    k = 0
    aliInfo = {}
    for k, fname in enumerate(fnames):
        prot = fname.split('.')[0]
        print prot
        infile = MEME_OUT_DIR+prot+'/meme.txt'
        if not os.path.exists(infile):
            continue
        info = parseMemeTxt(infile)
        mInfo = info['motifInfo'][0]
        #print len(mInfo['aliSeqs'])
        mInfo_rw = makeWeightedMotif(mInfo['aliSeqs'])
        mInfo_rw_log2 = makeWeightedMotif(mInfo['aliSeqs'], log2 = True)
        msout.write('\t'.join([prot]+[mInfo[x] \
                             for x in mSummLabs[1:]])+'\n')
        for fhand in [mout, mout_rw, mout_rw_log2]:
            fhand.write('TF\t%s\nTF Name\t%s\nMotif\tM%d_chu\n' %(prot,prot, k+1))
            fhand.write('\t'.join(['Pos','A','C','G','T']) + '\n')
        for j in range(len(mInfo['pfm'])):
            mout.write('\t'.join([str(j+1)]+mInfo['pfm'][j])+'\n')
            mout_rw.write('\t'.join([str(j+1)]+['%f\t%f\t%f\t%f\n'%tuple(mInfo_rw[j])]))
            mout_rw_log2.write('\t'.join([str(j+1)]+['%f\t%f\t%f\t%f\n'%tuple(mInfo_rw_log2[j])]))
        for fhand in [mout, mout_rw, mout_rw_log2]:
            fhand.write('\n\n')
        bsout.write('\t'.join([prot]+[info['letterFreq'][x] \
                                 for x in bSummLabs[1:-1]])+'\tletter\n')
        bsout.write('\t'.join([prot]+[info['bgFreq'][x] \
                                 for x in bSummLabs[1:-1]])+'\tbg\n')
    msout.close()
    bsout.close()
    mout.close()
    mout_rw.close()
    mout_rw_log2.close()

if __name__ == '__main__':
    main()
