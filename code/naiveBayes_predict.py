# A script for making predictions based on parameter
# estimates output by gibbAlign_naiveBayes.py

import numpy as np
import os, sys
from gibbsAlign_naiveBayes import getHomeoboxData, makeAllLogos
from gibbsAlign_naiveBayes import getAlignedPWMs, getOrientedPWMs, reverseBaseOrient
from gibbsAlign_naiveBayes import getCondModel, getLogLikelihood, getLLsum, sampleStartPos
from getHomeoboxConstructs import subsetDict
from scipy.stats import pearsonr
from copy import deepcopy
from matAlignLib import comp_matrices, matrix_compl, PCC

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



if len(sys.argv) > 1:
    TST_SET, CORE_POS, OBS_GRPS, APOS_CUT, EDGE_CUT = sys.argv[1:6]
    MWID, N_CHAINS, RAND_SEED = \
        [int(x) for x in sys.argv[7:]]
    MAX_EDGES_PER_BASE = sys.argv[6]
    if MAX_EDGES_PER_BASE == 'None':
        MAX_EDGES_PER_BASE = None
    else:
        MAX_EDGES_PER_BASE = int(MAX_EDGES_PER_BASE)
else:
    TST_SET = 'b08'
    CORE_POS = 'useStructInfo' #'canon9'#
    OBS_GRPS = 'grpIDcore' #'hd1' #'none' #'clusters' #
    APOS_CUT = 'cutAApos_1.0_0.05'#'cutAApos_X'#'cutAApos_0.0_0.05' #
    EDGE_CUT = 'edgeCut_1.0_0.05'#'cutEdge_X' #
    MAX_EDGES_PER_BASE = None #2 # None means to ignore paramter
    MWID = 6
    N_CHAINS = 100
    RAND_SEED = 382738372#78374223 #

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

def readParams_gibbs(fname_cond, fname_marg, logProbs = True):
    """ Read in parameter estimates output by gibbAlign_naiveBayes.py
    as logProbs.
    """

    fin = open(fname_cond, 'r')
    h = {x:i for i,x in enumerate(fin.readline().strip().split())}
    cond = {}
    for line in fin:
        l = line.strip().split('\t')
        ap, bp, aa, b, p = \
            int(l[h['apos']]), int(l[h['bpos']]), l[h['aa']], l[h['base']],\
            float(l[h['prob']])
        if not cond.has_key((ap,bp)):
            cond[(ap,bp)] = np.zeros((len(AMINO),len(BASE)),dtype = 'float')
        if logProbs:
            cond[(ap,bp)][A2IND[aa],B2IND[b]] = np.log(p)
        else:
            cond[(ap,bp)][A2IND[aa],B2IND[b]] = p
    fin.close()

    fin = open(fname_marg, 'r')
    h = {x:i for i,x in enumerate(fin.readline().strip().split())}
    marg = {}
    for line in fin:
        l = line.strip().split('\t')
        bp, b, p = \
            int(l[h['bpos']]), l[h['base']], float(l[h['prob']])
        if not marg.has_key(bp):
            marg[bp] = np.zeros(len(BASE), dtype = 'float')
        if logProbs:
            marg[bp][B2IND[b]] = np.log(p)
        else:
            marg[bp][B2IND[b]] = p
    fin.close()

    return cond, marg

def getOffsets(runCmpFile, optRun = 'maxll'):
    """ Returns the correct pwm offsets and orientations based
    on the cmpRuns.txt file output by gibbAlign_naiveBayes.py
    """

    # Find the run number with the optimal likelihood
    if optRun == 'maxll':
        fin = open(runCmpFile, 'r')
        h = {x:i for i,x in enumerate(fin.readline().strip().split())}
        maxRun = 0; maxll = -1e30
        for line in fin:
            ll = float(line.strip().split()[h['loglik']])
            if ll > maxll:
                maxll, maxRun = ll, int(line.strip().split()[h['nRun']])
        fin.close()
    else:
        maxRun = optRun

    # Read in the appropriate start and offsets and orientations
    fin = open(runCmpFile, 'r')
    h = {x:i for i,x in enumerate(fin.readline().strip().split())}
    start = {}
    rev = {}
    for line in fin:
        l = line.strip().split()
        if int(l[h['nRun']]) == maxRun:
            start[l[h['prot']]] = int(l[h['start']])
            rev[l[h['prot']]] = int(l[h['rev']])
            mWid = int(l[h['mWid']])
            reorient = int(l[h['reorient']])
    fin.close()

    return start, rev, mWid, reorient

def makeNBpreds(cond, marg, edges, featVect, f2i):
    """ Make a naive bayes pred based on the conditional matices
    and marginal probabilities (log space probs are expected).
    """

    preds = []
    for v in featVect:
        pred = np.zeros((len(marg), len(marg[0])), dtype = 'float')
        for i in range(len(pred)):
            for j, f in enumerate(v):
                if j in edges[i]:
                    pred[i,:] += cond[(j,i)][f2i[f],:]
            pred[i,:] += marg[i]
            pred[i,:] = np.exp(pred[i,:])
            pred[i,:] /= pred[i,:].sum()
        preds.append(pred)
    return preds

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
    # for the best IC-PCC alignment to the fly pwm?

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

def evaluateFit(pwms, core, edges, paramDir,
                optRun = 'maxll', makeLogos = False):
    """ Evaluate the fit of the model to the training data
    """

    # Read in the estimated parameters and make predictions
    cond, marg = readParams_gibbs(paramDir+'condProbs.txt',
                                  paramDir+'baseProbs.txt')
    preds = makeNBpreds(cond, marg, edges,
                        [core[k] for k in sorted(core.keys())], A2IND)
    preds = {k: preds[i] for i,k in enumerate(sorted(core.keys()))}

    if makeLogos:
        logoDir = paramDir + 'predLogos/'
        if not os.path.exists(logoDir):
            os.makedirs(logoDir)
        makeAllLogos(preds, core, logoDir)

    # Get the aligned and oriented experimental PWMs
    start, rev, mWid, reorient = getOffsets(paramDir + \
                                            'compareRuns/cmpRuns.txt',
                                            optRun = optRun)
    aliPWMS = getAlignedPWMs(getOrientedPWMs(pwms, rev),
                             core, start, mWid, flipAli = reorient)

    # Output per-column agreement information
    fitDir = paramDir + 'fitInfo/'
    if not os.path.exists(fitDir):
            os.makedirs(fitDir)
    fout = open(fitDir+'perColumnAlignedVsPred.txt', 'w')
    fout.write('prot\tbpos\tic.y\tic.y.hat\tpcc\trmse\n')
    for k in sorted(preds.keys()):
        for i in range(len(preds[preds.keys()[0]])):
            x, y = aliPWMS[k][i], preds[k][i]
            pcc, sqErr, icy, icy_hat = \
                pearsonr(x, y)[0], np.sqrt(((x-y)**2).mean()), ICbits(x), ICbits(y)
            #print pcc, sqErr, icy, icy_hat
            fout.write('%s\t%d\t%e\t%e\t%e\t%e\n' \
                       %(k, i, icy, icy_hat, pcc, sqErr))
    fout.close()

    return aliPWMS

def makeSimilarityMats(aliPWMS, core, paramDir):

    # Create similarity matrices between protein core sequences
    # and the alignment inferred by the model
    prots = sorted(core.keys())
    protSim = np.zeros((len(prots),len(prots)), dtype = 'float')
    aliSim = np.zeros((len(prots),len(prots)), dtype = 'float')
    for i, k1 in enumerate(prots):
        c1 = core[k1]
        p1 = aliPWMS[k1]
        for j, k2 in enumerate(prots):
            c2 = core[k2]
            p2 = aliPWMS[k2]
            hdist = 0
            for xpos in range(len(c1)):
                if c1[xpos] != c2[xpos]:
                    hdist += 1
            protSim[i,j] = len(c2) - hdist
            aliSim[i,j] = \
                np.mean([np.sqrt(((p1[k]-p2[k])**2).mean()) \
                        for k in range(len(p1))])
    #protSim = len(prots[0]) - protSim
    aliSim = 1 - aliSim/np.max(aliSim)

    fitDir = paramDir + 'fitInfo/'
    if not os.path.exists(fitDir):
            os.makedirs(fitDir)
    np.savetxt(fitDir+'coreSimilarities_hammDist.txt', protSim)
    np.savetxt(fitDir+'alignmentSimilarities_rmse.txt', aliSim)
    fout = open(fitDir+'protnames.txt', 'w')
    fout.write("\n".join(prots))
    fout.close()

    return protSim, aliSim

def main():
    tstSet = TST_SET

    # Create the to level output directory
    suffix = ''
    if OBS_GRPS != 'none':
        suffix += '_'+OBS_GRPS
    suffix += '/%s_rSeed_%d/mWid_%d/nChains_%d/' \
        %(tstSet, RAND_SEED, MWID, N_CHAINS)
    if CORE_POS == 'useStructInfo':
        suffix += APOS_CUT+'_'+EDGE_CUT
    elif CORE_POS == 'canon9':
        suffix += 'canon9_'+EDGE_CUT
    if MAX_EDGES_PER_BASE is None:
        suffix += '/'
    else:
        suffix += '_maxEdges%d/' %MAX_EDGES_PER_BASE
    paramDir = '../homeodomain/gibbsAlign_output%s' %(suffix)

    # Read in and filter the data the same way as in gibbsAlign_naiveBayes.py
    if CORE_POS == 'canon9':
        aaPosList = CANON9
    else:
        aaPosList = CORE_POS
    seqs, pwms, core, full, trunc, aaPosList, edges, edges_hmmPos = \
        getHomeoboxData(['n08', 'b08'], MWID, aaPosList = aaPosList,
                        aposCut = APOS_CUT, edgeCut = EDGE_CUT,
                        maxEdgesPerBase = MAX_EDGES_PER_BASE,
                        N51A_bpos = (MWID-6)/2+2)

    # For mapping/validation
    n08seqs, n08pwms, n08core = seqs['n08'], pwms['n08'], core['n08']

    # Subset to set of interest
    seqs, pwms, core, full, trunc = seqs[tstSet], pwms[tstSet], \
        core[tstSet], full[tstSet], trunc[tstSet]
    keep = set([k for k in core.keys() if 'X' not in core[k]])
    [subsetDict(x,keep) for x in [seqs, pwms, core, full, trunc]]

    print(paramDir)
    # Evaluate the fit of the aligned pwms to the inferred model
    aliPWMS = evaluateFit(pwms, core, edges, paramDir, optRun = 'maxll',
                          makeLogos = False)

    # Compute similary matrices for clustering test set
    protSim, aliSim = makeSimilarityMats(aliPWMS, core, paramDir)

    # Compare distance to nearest core/prot in Noyes 08 set
    prot2grp, grp2prot = getNoyesSpecGrps(NOYES_SPEC_GRPS)    
    closestCore = getClosestNoyesProt(core, n08core)
    aliScoresClosest = bestAliScoreWithNoyes(closestCore, aliPWMS, n08pwms)

    # Group aligned PWMs that have HD-0 cores with members of Noyes 08 set
    oldPath = paramDir+'logos_1/'
    newPath = paramDir+'logos_grpByNoyes08_hd0/'

    if not os.path.exists(newPath):
        os.makedirs(newPath)
    fout = open(paramDir+'n08_closestCore.txt', 'w')
    fout.write('prot\tprot.n08\tcore\tcore.n08\tgrp.n08\thDist\tbestAli\n')
    for k in sorted(closestCore.keys()):
        p, p_n08, dist, c, c_n08, g, aliScore= \
            k, closestCore[k][0], closestCore[k][1], core[k], \
            n08core[closestCore[k][0]], prot2grp[closestCore[k][0]], \
            aliScoresClosest[k]
        fout.write('%s\t%s\t%s\t%s\t%s\t%d\t%f\n' \
                   %(p, p_n08, c, c_n08, g, dist, aliScore))
        if dist == 0:
            oldLab, newLab = '_'.join([p,c]), '_'.join([g,p_n08,p,c])
            os.system('cp %s %s' %(oldPath+oldLab+'.pdf',
                                   newPath+newLab+'.pdf'))
    fout.close()

if __name__ == '__main__':
    main()
