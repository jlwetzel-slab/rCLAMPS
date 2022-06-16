# A script for making predictions based on parameter
# estimates output by gibbAlign_GLM.py [predict fly proteins using the model built using mouse proteins]

import numpy as np
import re
import os, sys
from gibbsAlign_GLM import getHomeoboxData, makeAllLogos
from gibbsAlign_GLM import getTestProtsAndPWMs, getTrainPairsAndInfo
from gibbsAlign_naiveBayes import getAlignedPWMs, getOrientedPWMs, reverseBaseOrient
from gibbsAlign_naiveBayes import getCondModel, getLogLikelihood, getLLsum, sampleStartPos
from getHomeoboxConstructs import subsetDict
from scipy.stats import pearsonr
from copy import deepcopy
from matAlignLib import comp_matrices, matrix_compl, PCC
import pickle

NOYES_SPEC_GRPS = '../homeodomain/flyFactorSurvey/noyes_cell08/'+\
    'specificity_groups.txt'
BASE = ['A','C','G','T']
REV_COMPL = {'A':'T','C':'G','G':'C','T':'A'}
AMINO = ['A','C','D','E','F','G','H','I','K','L',
         'M','N','P','Q','R','S','T','V','W','Y']
B2IND = {x: i for i, x in enumerate(BASE)}
A2IND = {x: i for i, x in enumerate(AMINO)}
IND2B = {i: x for i, x in enumerate(BASE)}
IND2A = {i: x for i, x in enumerate(AMINO)}

CANON9 = [2,3,5,6,47,50,51,54,55]
MAX_EDGES_PER_BASE = None
OBS_GRPS = 'grpIDcore' #'hd1' #'none' #'clusters' #
CORE_POS = 'useStructInfo'    # Uses the precomputed structural alignments
OBS_GRPS = 'grpIDcore'        # Perform group updates based on common "core" AAs in proteins
APOS_CUT = 'cutAApos_1.0_0.05'  # Controls AA contact threshold for binarization
EDGE_CUT = 'edgeCut_1.0_0.05'   # Controls AA contact threshold for binarization
TR_SET = 'cisbp'#'b08'  # Just controls where code looks for the PCMs

MWID = 6  # Number of base positions in the PDSIM
RAND_SEED = 382738375 #78374223 # Random seed for reproducibilty
N_CHAINS = 100

MODEL_FILE = '../results/cisbp-chuAll/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50.pickle'

def entropyBits(f):
    return -(f*np.log2(f)).sum()

def ICbits(f):
    return np.log2(len(f)) - entropyBits(f)

def getNoyesSpecGrps(fname):
    """ Read file return mapping between fly prot and specificity group
    as defined by clustering in Christensen/Noyes Cell 08
    """
    fin = open(fname, 'r')
    prot2grp = {}
    grp2prot = {}
    for line in fin:
        grp, prot = [x for x in line.strip().split('\t')]
        if not grp2prot.has_key(grp):
            grp2prot[grp] = [prot]
        else:
            grp2prot[grp].append(prot)
        prot2grp[prot] = grp
    fin.close()
    return prot2grp, grp2prot

def getClosestNoyesProt(seq, ref):
    """ Returns a dictionary mapping labels from dict seq to
    a (label, distance) pair corresponding to an entry from
    the ref dict that is closest in terms of hamming distance
    """
    closest = {}
    for k in seq.keys():
        minDist = 1e9
        minKey = ''
        x = seq[k]
        for r in ref.keys():
            dist = 0
            y = ref[r]
            for i in range(len(x)):
                if x[i] != y[i]:
                    dist += 1
            if dist < minDist:
                minKey = r
                minDist = dist
        closest[k] = (minKey, minDist)
    return closest

def bestAliScoreWithNoyes(closestCore, aliPWMS, n08pwms):
    # What is the average per-column pearson correlation
    # for the best alignment to the fly pwm?

    pccVal = {}
    for k1 in sorted(closestCore.keys()):#[:1]:
        k2, dist = closestCore[k1]
        m1, m2 = aliPWMS[k1], n08pwms[k2]
        score, shift, rev = comp_matrices(m2,m1,oneSided = True,
                                          minWidthM2 = len(m1))
        if rev:
            m2 = matrix_compl(m2)
        #print k1, k2, dist, score, shift, rev,
        i, j = 0, 0
        if shift > 0:
            j = shift
        elif shift < 0:
            i = -shift
        pccSum = 0.0
        while i<len(m1) and j<len(m2):
            pccSum += PCC(m1[i], m2[j])
            i, j = i+1, j+1
        pccVal[k1] = pccSum
        #print pccSum
    return pccVal

def main():

    with open(MODEL_FILE) as f:
        res = pickle.load(f)


    score = [x['ll'] for x in res]
    reorient = [x['reorient'] for x in res]
    start = [x['start'] for x in res]
    rev = [x['rev'] for x in res]
    opt = np.argmax(score)
    print(opt)

    trSet = TR_SET

    # Create the to level output directory
    suffix = ''
    if OBS_GRPS != 'none':
        suffix += '_'+OBS_GRPS
    suffix += '/cisbp-chu-scaled50/structFixed1_grpHoldout_glm_multinomial/%s_rSeed_%d/mWid_%d/nChains_%d/' \
        %(trSet, RAND_SEED, MWID, N_CHAINS)
    if CORE_POS == 'useStructInfo':
        suffix += APOS_CUT+'_'+EDGE_CUT
    elif CORE_POS == 'canon9':
        suffix += 'canon9_'+EDGE_CUT
    if MAX_EDGES_PER_BASE is None:
        suffix += '/'
    else:
        suffix += '_maxEdges%d/' %MAX_EDGES_PER_BASE
    paramDir = '../homeodomain/gibbsAlign_output%s' %(suffix)
    if not os.path.exists(paramDir):
            os.makedirs(paramDir)

    # Read in and filter the data the same way as in gibbsAlign_naiveBayes.py
    if CORE_POS == 'canon9':
        aaPosList = CANON9
    else:
        aaPosList = CORE_POS
    seqs, pwms, core, full, trunc, aaPosList, edges, edges_hmmPos = \
        getHomeoboxData(['n08', 'cisbp', 'chu'], MWID, aaPosList = aaPosList,
                        aposCut = APOS_CUT, edgeCut = EDGE_CUT,
                        maxEdgesPerBase = MAX_EDGES_PER_BASE,
                        N51A_bpos = (MWID-6)/2+2)

    # For mapping/validation
    n08pwms, n08core = pwms['n08'], core['n08']
    #cisseqs, cispwms, ciscore = seqs['cisbp'], pwms['cisbp'], core['cisbp']

    # Subset to set of interest
    pwms1, core1 = pwms[trSet], core[trSet]
    pwms2, core2 = pwms['chu'], core['chu']

    # Remove examples from cis-bp and chu that correspond to 
    # core seqs in the test set
    fin = open('../hd1-explore/0_splitChu2012-trainTest/testProts.txt','r')
    testProts = [x.strip() for x in fin.readlines()]
    fin.close()
    for x in [core1, core2, pwms1, pwms2]:
        allProts = x.keys()
        for prot in allProts:
            if prot in testProts:
                del x[prot]
    #print(len(pwms1), len(pwms2))

    # Combine pwms and core seqs remaining into a single dataset
    pwms, core = {}, {}
    for x in [pwms1, pwms2]:
        for prot in x.keys():
            pwms[prot] = x[prot]
    for x in [core1, core2]:
        for prot in x.keys():
            core[prot] = x[prot]
    #print len(pwms), len(core)

    keep = set([k for k in core.keys() if 'X' not in core[k]])
    [subsetDict(x,keep) for x in [pwms, core]]
    keep = set([k for k in core.keys() if '-' not in core[k]])
    [subsetDict(x, keep) for x in [pwms, core]]
    keep = set([k for k in pwms.keys() if pwms[k].shape[0] >= MWID])
    [subsetDict(x, keep) for x in [pwms, core]]

    print(paramDir)
    # Evaluate the fit of the aligned pwms to the inferred model
    flipAli = False
    if reorient[opt]:
        flipAli = True
    aliPWMS = getAlignedPWMs(getOrientedPWMs(pwms, rev[opt]), core, start[opt], 
                             MWID, flipAli = flipAli)
    logoDir = paramDir + 'logos_1/'
    if not os.path.exists(logoDir):
        os.makedirs(logoDir)
    makeAllLogos(aliPWMS, core, logoDir)


    # Compare distance to nearest core/prot in Noyes 08 set
    prot2grp, grp2prot = getNoyesSpecGrps(NOYES_SPEC_GRPS)
    closestCore = getClosestNoyesProt(core, n08core)
    aliScoresClosest = bestAliScoreWithNoyes(closestCore, aliPWMS, n08pwms)

    # Group aligned PWMs that have HD-0 cores with members of Noyes 08 set
    oldPath = paramDir+'logos_1/'
    newPath = paramDir+'logos_grpByNoyes08_hd0/'

    #print closestCore
    if not os.path.exists(newPath):
        os.makedirs(newPath)
    fnames = os.listdir(oldPath)
    fout = open(paramDir+'n08_closestCore.txt', 'w')
    fout.write('prot\tprot.n08\tcore\tcore.n08\tgrp.n08\thDist\tbestAli\n')
    for k in sorted(closestCore.keys()):
        try:
            p, p_n08, dist, c, c_n08, g, aliScore= \
                k, closestCore[k][0], closestCore[k][1], core[k], \
                n08core[closestCore[k][0]], prot2grp[closestCore[k][0]], \
                aliScoresClosest[k]
        except KeyError:
            print "Key not found"
            continue
        fout.write('%s\t%s\t%s\t%s\t%s\t%d\t%f\n' \
                   %(p, p_n08, c, c_n08, g, dist, aliScore))
        if dist == 0:
            oldLab = '_'.join([p,c])
            oldFile = [x for x in fnames if re.match(oldLab,x) is not None][0]
            newFile = '_'.join([g,p_n08])+'_'+oldFile
            os.system('cp %s %s' %(oldPath+oldFile,newPath+newFile))
    fout.close()

if __name__ == '__main__':
    main()
