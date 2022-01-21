# Helper functions

from runhmmer import *
from matAlignLib import *
from pwm import makeNucMatFile, makeLogo
from copy import deepcopy
import numpy as np
import os
from scipy import stats


NOYES08_INFO = '../homeodomain/flyFactorSurvey/noyes_cell08/' + \
    '1-s2.0-S009286740800682X-mmc2_edited.csv'
NOYES08_PWM_SRC = '../pwms/homeodomain/flyFactorSurvey_20180126.txt'
BERGER08_FASTA = '../homeodomain/uniprobe_homeodomains_mouse.fa'
BERGER08_PWM_SRC = '../pwms/homeodomain/uniprobe_20180126_Cell08/'
HBOX_HMM = '../pfamHMMs/Homeobox.hmm'
HBOX_HMM_LEN = 57
HMM_NAME = 'Homeobox'
HMM_OFFSET = 2
CORE_POS = [x - HMM_OFFSET for x in [2,3,5,6,47,50,51,54,55]]
LOGO_DIR = '../pwms/homeodomain/logos/'
ALIGNMENT_FILE = '../pwms/homeodomain/n08_b08_alignments_normScores.txt'
ALIGNMENT_CLUST_FILE = \
    '../pwms/homeodomain/n08_b08_aliClusts_b08bestExemplar.txt'
    #'../../pwms/homeodomain/n08_b08_aliClusts_kmeans10.txt'


MAKE_LOGOS = True
MAKE_LOGOS_CLUST_PWM = False
MAKE_LOGOS_CLUST_PROT = False
COMPUTE_ALIGNMENTS = False

def mergeDicts(d1, d2, dsetNames = None):
    # Merges two dictionaries indexed by strings into one combined
    # dict ... if there is a shared key across the dicts, then it is
    # turned into two keys with either '_1' or '_2' appended to the key
    # Also returns a map from each name to its original dsetName
    # according to the ordered list dsetNames (i,e., [d1name, d2name])

    d = {}
    dset = {}
    sharedKeys = set(d1.keys()) & set(d2.keys())
    for k in d1.keys():
        tmp = k
        if k in sharedKeys:
            k = k+'_1'
        d[k] = d1[tmp]
        if dsetNames is not None:
            dset[k] = dsetNames[0]
    for k in d2.keys():
        tmp = k
        if k in sharedKeys:
            k = k+'_2'
        d[k] = d2[tmp]
        if dsetNames is not None:
            dset[k] = dsetNames[1]
    return d, dset

def subsetDict(d, subset = set()):
    for k in d.keys():
        if k not in subset:
            del d[k]

def readUniprobePWM(fname, lskip = 3):
    # Return a pwm from a uniprobe-style file

    fin = open(fname, 'r')
    [fin.readline() for i in range(lskip)]
    p = []
    for line in fin:
        p.append([float(x) for x in line.strip().split(':')[1].split()])
    p = norm_matrix(np.array(p).T)
    return p

def getUniprobePWMs(fdir, subset = set()):
    # Parses a directory of uniprobe pwm files file to get the subset of
    # PWMs corresponding to proteins listed in subset.
    # whichPWM corresponds to which pwm version we're interested
    # in or which PBM processing version (i.e., '.bml.' ...).
    # Returns a dictionary mapping protein names to pwms

    fnames = os.popen('ls '+fdir, 'r')
    pwms = {}
    for fname in fnames:
        fname = fname.strip()
        if fname == 'readme.txt' or '_secondary' in fname \
            or '_RC' in fname or '.bml' in fname:
            continue
        prot = fname.split('_')[0]
        if not pwms.has_key(prot):
            pwms[prot] = readUniprobePWM(fdir+fname, lskip = 3)
    return pwms

def parseNoyes08Table(fpath):
    # Returns a mapping from construct (or gene) name to
    # the aa sequence used by Noyes 2008 Cell B1H homeobox
    # article

    fin = open(fpath, 'r')
    aaNames = fin.readline().strip().split(',')
    fin.readline()
    aaSeqs = fin.readline().strip().split(',')
    return {aaNames[i]: aaSeqs[i] for i in range(len(aaNames))}

def getFlyFactorPWMs(fpath, subset = set(),
                     whichPWM = 'SOLEXA', countSmooth = 0):
    # Parses a flyFactor (Type B) flat file to get the subset of
    # PWMs corresponding to proteins listed in subset.
    # whichPWM corresponds to which pwm version we're interested
    # in (i.e., Cell, SOLEXA, etc ...).
    # Returns a dictionary mapping genes names to pwms

    pwms = {}
    fin = open(fpath, 'r')
    line = fin.readline()
    prevLabel = ''
    prevGene = ''
    i = 1
    while line != '' and line != '\n':
        if line[0] == '>':
            label = line.strip()
            l = line[1:-1].split('_')
            gene, pwmType = l[0], l[1]
            if pwmType != whichPWM:
                line = fin.readline()
                while line[0] != '>':
                    line = fin.readline()
                    i+=1
                    continue
            if prevLabel != '':
                pwms[prevGene] = norm_matrix(np.array(p))
            prevLabel = label
            prevGene = gene
            p = []
        else:
            p.append([float(x)+countSmooth \
                     for x in line.strip().split()])
        line = fin.readline()
        i+=1
    fin.close()
    pwms[prevGene] = norm_matrix(np.array(p))
    return pwms

def readFromFasta(fpath):
    # Reads seqs in a fasta file to a dict
    fin = open(fpath)
    line = fin.readline()
    prevLab = ''
    seqs = {}
    while line != '' and line != '\n':
        if line[0] == '>':
            if prevLab != '':
                seqs[prevLab] = s
            lab = line[1:-1]
            prevLab = lab
            s = ''
        else:
            s += line.strip()
        line = fin.readline()
    fin.close()
    seqs[lab] = s
    return seqs

def writeToFasta(seqDict, fpath):
    fout = open(fpath, 'w')
    for k in sorted(seqDict.keys()):
        fout.write('>%s\n' %k)
        fout.write('%s\n' %seqDict[k])
    fout.close()

def makeMatchStateTab(inpath, outpath, origSeqs, hmmLen, hmmName,
                      corePos = None):
    # Takes as input the output file created by parsehmmer*() and
    # convert to a table where for each match, each match state
    # has its own column, and its value is the amino acid occupying
    # that match state.
    # Returns a dictionary mapping strings to core sequences
    # according to the set of indices given by corePos

    # Default to using all aa positions
    if corePos is None:
        corePos = set(range(hmmLen))

    fin = open(inpath, 'r')
    seqs = {} # Full extended match states for domain HMM
    coreSeqs = {} # Concatenated "core" of extended match states for domain HMM
    truncated = {} # True if match was truncated, False otherwise
    scores = {}    # For book-keeping so we only keep the single best hit
    fin.readline() # Skip the header
    outStr = {}    # One line per protein
    for line in fin:
        l = line.strip().split()
        prot, eValue, bitScore, hmmStart, hmmEnd, protStart, \
            protEnd, ali, m = \
            l[0], l[2], float(l[3]), int(l[4]), int(l[5]), int(l[6]), int(l[7]), \
            l[8], l[9]

        if prot in scores.keys() and scores[prot] >= bitScore:
            continue

        alilen, mlen = len(ali), len(m)
        assert(alilen == mlen)

        # Map to match states
        i = 0  # index to current match sequence position
        j = 1  # index of the current match state in hmm
        truncated[prot] = False
        while i < alilen or j <= hmmLen:
            #print prot
            if i == 0 and j == 1:  # beginnign of line is protein name
                outStr[prot] = prot
            if j < hmmStart:  # HMM match is front truncated
                truncated[prot] = True
                protInd = protStart-(hmmStart-j+1)
                if protInd < 0:
                    outStr[prot] += ' X' # Truncation due to incomplete construct
                else:
                    outStr[prot] += ' %s' %origSeqs[prot][protInd]
                j += 1
            elif i >= len(ali): # HMM match is end truncated
                truncated[prot] = True
                protInd = protEnd-(hmmEnd-j+1)
                if protInd >= len(origSeqs[prot]):
                    outStr[prot] += ' X' # Truncation due to incomplete construct
                else:
                    outStr[prot] += ' %s' %origSeqs[prot][protInd]
                j += 1
            elif ali[i] == '.':  # insetion between match states
                i += 1
            elif m[i] == '.': # skipped match state
                outStr[prot] += ' X'
                i += 1
                j += 1
            else:       # Standard match state
                outStr[prot] += ' %s' %m[i]
                i += 1
                j += 1
        outStr[prot] += '\n'
        assert (len(outStr[prot]) - len(prot) - 1) == 2 * hmmLen
        seqs[prot] = ''.join([x for i, x in \
                             enumerate(outStr[prot].strip().split()[1:])])
        coreSeqs[prot] = \
            ''.join([x for i, x in
                enumerate(outStr[prot].strip().split()[1:]) if i in corePos])
        scores[prot] = bitScore
    fin.close()

    fout = open(outpath, 'w')
    for prot in sorted(outStr):
        fout.write(outStr[prot])
    fout.close()

    return coreSeqs, seqs, truncated

def getHDkNeighbors(seq, dbase, k = 1):
    # Given a sequence, seq, and a database of sequences as a
    # dict (dbase), returns the list of keys to items in dbase
    # that are hamming distance <= k from seq

    neighbors = {}
    for key in dbase.keys():
        s = dbase[key]
        #mm = 0   # number of mismatches
        if s == seq:
            neighbors[key] = []  # A hamming distance 0 match
            continue
        for i in range(len(seq)): # Check hamming dist
            if s[i] != seq[i]:
                #mm += 1
                if not neighbors.has_key(key):
                    neighbors[key] = []
                neighbors[key].append((i, s[i]))
                if len(neighbors[key]) > k:
                    del (neighbors[key])
                    break
    return neighbors

def getAllHDk_dbase(dbase1, dbase2, k = 1):
    # Given two databases mapping keys to sequences
    # returns a nested dictionary mapping keys from
    # dbase1 to keys in dbase2 for which these keys
    # map to seqs at most hamming distance k from
    # the seq keyed by the element in dbase1

    neighbors = {}
    for key in dbase1.keys():
        s = dbase1[key]
        neighbors[key] = getHDkNeighbors(s, dbase2, k = k)
    return neighbors

def makeNbrTable(nbrs, nbr_map, outpath, maxK = 1):

    fout = open(outpath, 'w')
    fout.write("id1\tseq1\tid2\tsubPos\tsubVal\thammDist\n")
    for k1 in nbrs.keys():
        seq1 = nbr_map[k1]
        for k2 in nbrs[k1].keys():
            if k1 == k2:
                continue
            subList = nbrs[k1][k2]
            if len(subList) > maxK:
                continue
            if len(subList) == 0:   # Handle the excat match case
                fout.write('%s\t%s\t%s\t%s\t%s\t%d\n' \
                           %(k1,seq1,k2,'NA','NA',len(subList)))
            for sub in subList:
                subPos, subVal = sub[0], sub[1]
                fout.write('%s\t%s\t%s\t%d\t%s\t%d\n' \
                           %(k1,seq1,k2,subPos,subVal,len(subList)))
    fout.close()

def makeAllLogosDset(pwms, coreSeqs, dset, logoDir, rc = False):
    # Place logos for every PWM in logoDir

    for gene in set(pwms.keys()) & set(coreSeqs.keys()):
        pwm, core = pwms[gene], coreSeqs[gene]
        if rc:
            pwm = matrix_compl(pwm)
        logoLab = '_'.join([gene, core, dset])
        makeNucMatFile('./tmp/','tmp',pwm)
        makeLogo('./tmp/tmp.txt',logoDir+logoLab+'.pdf',
                 alpha = 'dna', colScheme = 'classic')
        os.system('rm %s' %'./tmp/tmp.txt')

def makeAllLogosByClust(pwms, clustNums, coreSeqs,
                        dsets, logoDir):
    # Place logos for every PWM in logoDir

    for i in range(2):
        for gene in set(pwms.keys()) & set(coreSeqs.keys()):
            clust, pwm, core, dset = \
                int(clustNums[gene]), pwms[gene], coreSeqs[gene], dsets[gene]
            if i == 0:
                ld = logoDir + 'byClust/'
                logoLab = '_'.join(['%02d'%clust, core, gene, dset])
            elif i == 1:
                ld = logoDir + 'byCoreSeq/'
                logoLab = '_'.join([core, gene, dset])
            makeNucMatFile('./tmp/','tmp',pwm)
            makeLogo('./tmp/tmp.txt',ld+logoLab+'.pdf',
                     alpha = 'dna', colScheme = 'classic')
            os.system('rm %s' %'./tmp/tmp.txt')

def starAlignAllPWMs(pwms, anchor, anchorName):
    # Align all PWMs in the dataset pwms to the "anchor" pwm

    ali = {}
    npwms = {}
    for p in pwms.keys():
        pwm = pwms[p]
        if p == anchorName:
            npwms[p] = deepcopy(anchor)
            ali[p] = (4.0, 0, True)
            continue
        score, shift, rev = comp_matrices(anchor, pwm, 'PCC',
                                          normScore = True)
        ali[p] = (score, shift, rev)
        if rev:
            pwm = matrix_compl(pwm)
        npwm = np.zeros((len(anchor),4), dtype = 'float')
        j = -shift
        for i in range(len(anchor)):
            if j < 0 or j >= len(anchor) or j >= len(pwm):
                j += 1
                continue
            else:
                npwm[i] = pwm[j]
            j += 1
        npwms[p] = npwm
    return npwms, ali

def cmpSpecs(x, y):
    # Compares the columns of two 2D dicts of equal width

    pcc = {}
    for i in range(len(x)):
        pcc[i] = {}
        pcc[i]['pcc'] = stats.pearsonr(x[i], y[i])[0]
    return pcc


def cmpCoreSpecs(p1, c1, p2 = None, c2 = None):
    # Compare specificities for identical core sequences
    # within and across each group

    # Within set p1

    # Within set p2

    # Across the two sets
    pass

"""
    def compareSpecs(preds, specs):
    # Compares a set of predicted specificities to
    # a set of known specificities
    # - preds is a nested dictionary of predicted
    # per-base-position specs
    # - specs is the analogous dictionary of known specs
    # Returns a nested dictionary of per-base-position
    # measures of similarity between each base position
    # of the predicted and known specificity

    comp = {}
    for core in preds.keys():
        comp[core] = {}
        pred = preds[core]
        known = specs[core]
        for bpos in sorted(list(set(specs[core].keys()) \
                                & set(preds[core].keys()))):
            comp[core][bpos] = {}
            pred_col = pred[bpos]
            known_col = known[bpos]
            comp[core][bpos]['pcc'] = \
                stats.pearsonr(pred_col, known_col)[0]

    return comp

def outputComparisons(comp, fpath):
    # Outputs a table of comparisons that have
    # been computed by compareSpecs

    fout = open(fpath, 'w')

    # Make header
    headerCols = ['coreSeq', 'bpos']
    firstCore = comp.keys()[0]
    firstBpos = comp[firstCore].keys()[0]
    headerCols = headerCols + sorted(comp[firstCore][firstBpos].keys())
    fout.write(headerCols[0])
    for x in headerCols[1:]:
        fout.write('\t%s' %x)
    fout.write('\n')

    # Write the rows
    for core in sorted(comp.keys()):
        for bpos in sorted(comp[core].keys()):
            for measure in sorted(comp[core][bpos].keys()):
                score = comp[core][bpos][measure]
                fout.write("%s\t%d\t%f\n" %(core, bpos, score))

    fout.close()
"""

def alignAllPairsPWMs(pwms, id2seq, id2dset, writeTab = None,
                      normScore = False):
    # Returns a dictionary with the alignment information
    # for all pairs of pwms from the dataset pwms.

    ali = {}
    for i, p1 in enumerate(pwms.keys()):
        ali[p1] = {}
        for j, p2 in enumerate(pwms.keys()):
            print(i, j, p2)
            ali[p1][p2] = comp_matrices(pwms[p1], pwms[p2], 'PCC',
                                        normScore = normScore)

    if writeTab is not None:
        fout = open(writeTab, 'w')
        fout.write("id1\tseq1\tid2\tseq2\taliShift\taliRev\taliScore\t" +\
                   "dset1\tdset2\n")
        for p1 in ali.keys():
            for p2 in ali[p1].keys():
                score, shift, rev = ali[p1][p2]
                rev = int(rev)
                outStr = '%s\t%s\t%s\t%s\t%d\t%d\t%e\t%s\t%s\n' \
                    %(p1, id2seq[p1], p2, id2seq[p2], shift, rev, \
                      score, id2dset[p1], id2dset[p2])
                fout.write(outStr)
        fout.close()
    return ali

def getClusteredAliPWMs(pwms, clustAliTab, name2coreSeq, name2dset,
                        makeLogos = None):
    # Takes a set of pwms and precomputed alignment information
    # (based on clustering all-pairs alignment) and returns a
    # set of PWMs shifted according to this alignmnt information.
    # Optionally, creates logos according to these alignments
    # by supplying the argument with a folder name in which to
    # create the logos

    fin = open(clustAliTab, 'r')
    head = {x: i for i, x in enumerate(fin.readline().strip().split('\t'))}
    npwms = {}
    clustNums = {}
    for line in fin:
        l = line.strip().split()
        clust, hubProt, p, shift, rev = \
            int(l[head['clust.id1']]), l[head['id1']], l[head['id2']], \
            int(l[head['aliShift']]), int(l[head['aliRev']])
        hubPWM = pwms[hubProt]
        pwm = pwms[p]
        if p == hubProt:
            npwms[p] = deepcopy(hubPWM)
            clustNums[p] = clust
            continue
        if rev:
            pwm = matrix_compl(pwm)
        npwm = np.zeros((len(hubPWM),4), dtype = 'float')
        j = -shift
        for i in range(len(hubPWM)):
            if j < 0 or j >= len(hubPWM) or j >= len(pwm):
                j += 1
                continue
            else:
                npwm[i] = pwm[j]
            j += 1
        npwms[p] = npwm
        clustNums[p] = clust

    if makeLogos is not None:
        makeAllLogosByClust(npwms, clustNums, name2coreSeq,
                            name2dset, makeLogos)

    return npwms

def getClusterProtSeqLogos(clustAliTab, protSeqs, makeLogos = None):
    # Make logos based on protein sequences associated with
    # the clustered PWMs

    fin = open(clustAliTab, 'r')
    head = {x: i for i, x in enumerate(fin.readline().strip().split('\t'))}
    seqs = {}  # Lists of protein sequences keyed by cluster
    for line in fin:
        l = line.strip().split()
        clust, hubProt, p = \
            int(l[head['clust.id1']]), l[head['id1']], l[head['id2']]
        if not seqs.has_key((clust, hubProt)):
            seqs[(clust, hubProt)] = []
        seqs[(clust, hubProt)].append((p, protSeqs[p]))

    #print seqs
    #"""
    if makeLogos is not None:
        logoDir = makeLogos
        for clust, hubProt in seqs.keys():
            tmpFile = './tmp/tmp.fatsa'
            fout = open(tmpFile, 'w')
            numSeqs = len(seqs[(clust, hubProt)])
            for i, (p, s) in enumerate(seqs[(clust, hubProt)]):
                fout.write('>%s\n%s\n'%(p, s))
            fout.close()
            logoLab = '%02d_%s_nSeqs%03d'%(clust, hubProt, numSeqs)
            #print logoDir+logoLab+'.pdf'
            makeLogo(tmpFile,logoDir+logoLab+'.pdf',
                     datatype = 'fasta', alpha = 'protein',
                     colScheme = 'chemistry', size = 'large')
            os.system('rm %s' %tmpFile)
    #"""
    return seqs

def main():

    # Extract the seqs from the Cell08 paper and Berger08 paper
    noyes08Seqs = parseNoyes08Table(NOYES08_INFO)
    berger08Seqs = readFromFasta(BERGER08_FASTA)

    # Associate each homeodomain with its PWM
    n08_pwm = getFlyFactorPWMs(NOYES08_PWM_SRC,
                               subset = set(noyes08Seqs.keys()),
                               whichPWM = 'Cell', countSmooth = 1)
    b08_pwm = getUniprobePWMs(BERGER08_PWM_SRC,
                              subset = set(berger08Seqs.keys()))

    # Subset seqs to those with PWMs
    subsetDict(noyes08Seqs, set(n08_pwm.keys()))
    subsetDict(berger08Seqs, set(b08_pwm.keys()))

    # Create fastas to give to hmmer
    noyes08Fasta = '/'.join(NOYES08_INFO.split('/')[:-1]) + \
        '/ffs_homeodomains_fly_hasPWM.fa'
    writeToFasta(noyes08Seqs, noyes08Fasta)
    berger08Fasta = '/'.join(BERGER08_FASTA.split('/')[:-1]) + \
        '/uniprobe_homeodomains_mouse_hasPWM.fa'
    writeToFasta(berger08Seqs, berger08Fasta)

    # Run hmmer on fasta files
    noyes08Descs = getdescs(noyes08Fasta)
    berger08Descs = getdescs(berger08Fasta)
    noyes08HmmerOut = '/'.join(NOYES08_INFO.split('/')[:-1]) + \
        '/ffs_homeodomains_fly_hasPWM.hmmer3.out.txt'
    berger08HmmerOut = '/'.join(berger08Fasta.split('/')[:-1]) + \
        '/uniprobe_homeodomains_mouse_hasPWM.hmmer3.out.txt'
    runhmmer3(noyes08HmmerOut,HBOX_HMM,noyes08Fasta,HMM_NAME,noyes08Descs)
    runhmmer3(berger08HmmerOut,HBOX_HMM,berger08Fasta,HMM_NAME,berger08Descs)

    # Make match state table for each dataset and extract core seqs
    noyes08MatchTab = '/'.join(NOYES08_INFO.split('/')[:-1]) + \
        '/ffs_homeodomains_fly_hasPWM.matchTab.txt'
    berger08MatchTab = '/'.join(berger08Fasta.split('/')[:-1]) + \
        '/uniprobe_homeodomains_mouse_hasPWM.matchTab.txt'
    n08, n08_full, n08_trunc = \
        makeMatchStateTab(noyes08HmmerOut, noyes08MatchTab, noyes08Seqs,
                          HBOX_HMM_LEN, HMM_NAME, corePos = CORE_POS)
    b08, b08_full, b08_trunc = \
        makeMatchStateTab(berger08HmmerOut, berger08MatchTab,berger08Seqs,
                          HBOX_HMM_LEN, HMM_NAME, corePos = CORE_POS)

    # Subset pwms to those which have hmmer3 matches
    # **** (Shouldn't they *all* have matches though??
    #  Somehow one missing from each dataset ... come back to this later)
    subsetDict(n08_pwm, set(n08.keys()))
    subsetDict(b08_pwm, set(b08.keys()))

    # Get the hd-k neighbors and create tables
    k = 8
    n08_nbr = getAllHDk_dbase(n08, n08, k = k)
    b08_nbr = getAllHDk_dbase(b08, b08, k = k)
    n08_b08_nbr = getAllHDk_dbase(n08, b08, k = k)
    n08_nbr_out = '.'.join(noyes08MatchTab.split('.')[:-2])+'.nbrsHD%d.txt'%k
    b08_nbr_out = '.'.join(berger08MatchTab.split('.')[:-2])+'.nbrsHD%d.txt'%k
    n08_b08_nbr_out = '../homeodomain/n08_v_b08_hasPWM.nbrsHD8.txt'
    makeNbrTable(n08_nbr, n08, n08_nbr_out, maxK = k)
    makeNbrTable(b08_nbr, b08, b08_nbr_out, maxK = k)
    makeNbrTable(n08_b08_nbr, n08, n08_b08_nbr_out, maxK = k)

    # Prior to aligning, let's visually look at the PWMs as keyed
    # by their core sequences with original registers
    if MAKE_LOGOS:
        makeAllLogosDset(n08_pwm, n08, 'n08', LOGO_DIR+'noyes08/')
        makeAllLogosDset(n08_pwm, n08, 'n08', LOGO_DIR+'noyes08_rc/',
                         rc = True)
        makeAllLogosDset(b08_pwm, b08, 'b08', LOGO_DIR+'berger08/')
        makeAllLogosDset(b08_pwm, b08, 'b08', LOGO_DIR+'berger08_rc/',
                         rc = True)


    print(ALIGNMENT_FILE)
    # Combine the sets of PWMs and names to all-pairs alignments
    combinedPWMs, _ = mergeDicts(n08_pwm, b08_pwm)
    combinedCores, dsetMap = mergeDicts(n08, b08,
                                        dsetNames = ['n08','b08'])
    combinedFull, dsetMap = mergeDicts(n08_full, b08_full,
                                        dsetNames = ['n08','b08'])

    if COMPUTE_ALIGNMENTS:
        # Compute all pairwise alignments then cluster
        # based on alignment scores.
        aliInfo = alignAllPairsPWMs(combinedPWMs, combinedCores,
                                    dsetMap, writeTab = ALIGNMENT_FILE,
                                    normScore = True)
    else:
        clustDir = LOGO_DIR+'clusters/%s/'\
            %ALIGNMENT_CLUST_FILE[:-4].split('_')[-1]+'/'
        if MAKE_LOGOS_CLUST_PWM:
            pwmDir = clustDir
        else:
            pwmDir = None
        if MAKE_LOGOS_CLUST_PROT:
            fullDir = clustDir+'protSeqClustLogos/'
            coreDir = clustDir+'protSeqClustLogos_core/'
        else:
            fullDir = None
            coreDir = None
        aliPWM = getClusteredAliPWMs(combinedPWMs, ALIGNMENT_CLUST_FILE,
                                     combinedCores, dsetMap,
                                     makeLogos = pwmDir)
        seqsByClust = getClusterProtSeqLogos(ALIGNMENT_CLUST_FILE, combinedFull,
                                             makeLogos = fullDir)
        coresByClust = getClusterProtSeqLogos(ALIGNMENT_CLUST_FILE, combinedCores,
                                              makeLogos = coreDir)

if __name__ == '__main__':
    main()
