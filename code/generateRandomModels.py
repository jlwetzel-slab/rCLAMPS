# Generate random models to calculate their likelihoods

# Uses an alignemnt from gibbAlign_GLM.py to make predicitons in a 
# holdout validation scenario where for each held out protein,
# no other protein with and identical set base-contacting amino acid
# residues can be included in the training set.

import numpy as np
import scipy
import os, sys
from pwm import makeNucMatFile, makeLogo
from matAlignLib import comp_matrices, matrix_compl, information_content
from gibbsAlign_GLM import getHomeoboxData, makeAllLogos, makePWMtab
from gibbsAlign_GLM import getAlignedPWMs, getOrientedPWMs
from gibbsAlign_GLM import assignObsGrps, formGLM_fullX, formGLM_testX
from gibbsAlign_GLM import formGLM_trainX, formGLM_trainW, formGLM_Y, createGLMModel
from gibbsAlign_GLM import getTrainPairsAndInfo, initStarts, form_model, readSeedAlignment, computeGLMLoglikelihood
from getHomeoboxConstructs import subsetDict
import pickle
import time

# Used by various functions - do not change these
BASE = ['A','C','G','T']
REV_COMPL = {'A':'T','C':'G','G':'C','T':'A'}
AMINO = ['A','C','D','E','F','G','H','I','K','L',
         'M','N','P','Q','R','S','T','V','W','Y']
B2IND = {x: i for i, x in enumerate(BASE)}
A2IND = {x: i for i, x in enumerate(AMINO)}
IND2B = {i: x for i, x in enumerate(BASE)}
IND2A = {i: x for i, x in enumerate(AMINO)}
COMPL_BASE = {0:3, 1:2, 2:1, 3:0}

CORE_POS = 'useStructInfo'    # Uses the precomputed structural alignments
OBS_GRPS = 'grpIDcore'        # Perform group updates based on common "core" AAs in proteins
APOS_CUT = 'cutAApos_1.0_0.05'  # Controls AA contact threshold for binarization
EDGE_CUT = 'edgeCut_1.0_0.05'   # Controls AA contact threshold for binarization
MAX_EDGES_PER_BASE = None       # None means to just ignore this parameter
TRAIN_SET = 'cisbp'#'b08'
MWID = 6  # Number of base positions in the PDSIM
RESCALE_PWMS = True
EXCLUDE_TEST = False
MAKE_LOGOS = False
ADD_RESIDUES = False         # if True include additional residues from Christensen et al. (2012)
SEED_FILE = '../precomputedInputs/fixedStarts_homeodomains_all.txt'

def makePCCtable(exp, pred, core, fname):
    # Compute the per-position pearson correlation coefficients
    # and mean squared errors between the predicted and the 
    # experimental pwms assuming that the alignment s is correct

    fout = open(fname, 'w')
    fout.write('prot\tpos\tpcc\trmse\tic.exp\tic.pred\tcoreSeq\n')
    for prot in exp.keys():
        for pos in range(len(exp[prot])):
            x, y = exp[prot][pos], pred[prot][pos]
            pcc = scipy.stats.pearsonr(x,y)[0]
            rmse = np.sqrt(np.sum((x - y)**2))
            icPred, icExp = information_content(x), information_content(y)
            fout.write('%s\t%d\t%f\t%f\t%f\t%f\t%s\n'\
                       %(prot, pos+1, pcc, rmse, icPred, icExp, core[prot]))
    fout.close()

def main():

    if EXCLUDE_TEST:
        dirStem = '../results/cisbp-chu/'
        labStem = 'cisbp-chu/'
    else:
        dirStem = '../results/cisbp-chuAll/'
        labStem = 'cisbp-chuAll/'

    #mainOutDir = dirStem+'structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50'

    mainOutDir = '../my_results/allHomeodomainProts/'

    # Obtain Model
    #filename = dirStem+'structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50.pickle'
    filename = mainOutDir+'result.pickle'
    with open(filename) as f:
        res = pickle.load(f)
    score = [x['ll'] for x in res]
    opt = np.argmax(score)
    reorient = [x['reorient'] for x in res][opt]
    start = [x['start'] for x in res][opt]
    rev = [x['rev'] for x in res][opt]
    model = [x['final_model'] for x in res][opt]
    print('opt: %d' %opt)


    # Read data
    if CORE_POS == 'canon9':
        aaPosList = CANON9
    else:
        aaPosList = CORE_POS

    trainPWMs, trainCores, full, edges, edges_hmmPos, aaPosList, testProts = \
        getTrainPairsAndInfo(rescalePWMs = RESCALE_PWMS, excludeTestSet = EXCLUDE_TEST,
                             addResidues = ADD_RESIDUES)

    print edges
    print edges_hmmPos
    print aaPosList
    #print trainCores

    obsGrps = assignObsGrps(trainCores, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots  

    nDoms = {}
    for p in uniqueProteins:
        nDoms[p] = len(trainCores[p])/len(aaPosList)

    fixedStarts = readSeedAlignment(SEED_FILE, include = trainPWMs.keys())

    fout = open(mainOutDir+'randomModelScores.txt', 'w')
    fout.write('score\n')
    fullX, grpInd = formGLM_fullX(trainCores, edges, uniqueProteins, obsGrps)
    for k in range(100):
        start, rev = initStarts(uniqueProteins, trainPWMs, len(edges.keys()), nDoms, fixedStarts = fixedStarts)
        model = form_model(fullX, uniqueProteins, nDoms, trainPWMs, start, rev)
        ll = 0
        currInd = 0
        # Cycle through all the domains to compute total GLM log-likelihood
        for p in uniqueProteins:
            nextInd = currInd + 4*nDoms[p]
            testX = formGLM_testX(fullX, currInd, nextInd-1)
            testW = formGLM_trainW(trainPWMs, [p], nDoms, start, rev)
            for i, ind in enumerate(range(currInd, nextInd, 4)):
                testX_dom = formGLM_testX(fullX, ind, ind+4-1)
                testW_dom = {}
                for j in range(MWID):
                    # Each domain in a multidomain protein contributes equally to the likelihood function
                    testW_dom[j] = testW[j][i*4:(i+1)*4,]/sum(testW[j][i*4:(i+1)*4,])
                ll += computeGLMLoglikelihood(testX_dom, testW_dom, model)
            currInd = nextInd
        fout.write('%f\n'%ll)
        print k, ll
    fout.close()

if __name__ == '__main__':
    main()