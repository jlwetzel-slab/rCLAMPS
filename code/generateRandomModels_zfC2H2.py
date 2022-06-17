# Uses an alignemnt from gibbAlign_GLM.py to make predicitons in a 
# holdout validation scenario where for each held out protein,
# no other protein with and identical set base-contacting amino acid
# residues can be included in the training set.

import numpy as np
import scipy, os, sys
from pwm import makeNucMatFile, makeLogo
from matAlignLib import comp_matrices, matrix_compl, information_content
from gibbsAlign_GLM import getPrecomputedInputs_zfC2H2, makeAllLogos, makePWMtab, readSeedAlignment
from gibbsAlign_GLM import assignObsGrps, formGLM_fullX, formGLM_testX
from gibbsAlign_GLM import formGLM_trainX, formGLM_trainW, formGLM_Y, createGLMModel, initStarts, form_model, computeGLMLoglikelihood
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

MWID = 4  # Number of base positions in the PDSIM
RIGHT_OLAP = 1
MAKE_LOGOS = True
OBS_GRPS = 'grpIDcore'
ANCHOR_B1H = False

# Input files for PWMs and protein info
PROT_SEQ_FILE = '../precomputedInputs/zf-C2H2/prot_seq_fewZFs_hmmerOut_clusteredOnly_removeTooShort.txt'  # Input protein domain file subsetted to relvant amino acid contacting positions
PROT_SEQ_FILE_FFS = '../flyFactorSurvey/enuameh/enuameh_perFinger_processedProtInfo.txt'
PWM_INPUT_TABLE = '../precomputedInputs/zf-C2H2/pwmTab_fewZFs_clusteredOnly_removeTooShort.txt'   # A table of PWMs corresponding to prots in PROT_SEQ_FILE
PWM_INPUT_FILE_FFS = '../flyFactorSurvey/enuameh/flyfactor_dataset_A.txt'
SEED_FILE = '../flyFactorSurvey/enuameh/enuameh_startPosInfo.txt'
DIVERSE_ZF_SET = ['bowl','CG31670','ken','ovo','pho','Sp1']  # Set to None to seed with random ZFs instead

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

def predictSpecificity_array_ZF(fullX, model, startInd, arrayLen):
    # Predicts the specificity for an array of tandem domains with potential 
    # overlap in DBD-base interfaces (e.g., for zf-C2H2 proteins)

    # Predict an MWID length specificity for each domain individually
    testX = formGLM_testX(fullX, startInd, (startInd+4*arrayLen)-1)
    #print testX[0].shape
    pwms = {}
    for i in range(arrayLen):
        pwm = []
        for j in range(MWID):
            prediction = model[j].predict_proba(testX[j][i*4:(i+1)*4,])[0].tolist()
            pwm.append(prediction)
        pwms[i] = pwm
    
    #for k in sorted(pwms.keys()):
    #    print pwms[k]

    # Average the overlapping positions for adjacent domains and concatenate
    pwm = []
    for i in range(arrayLen):
        if i == 0:
            p = np.array(pwms[i])
            p[MWID-1,:] = (p[MWID-RIGHT_OLAP,:] + np.array(pwms[i+1][0]))/float(2)
            pwm += p.tolist()            
        elif i == arrayLen-1:
            pwm += pwms[i][RIGHT_OLAP:]
        else:
            p = np.array(pwms[i][RIGHT_OLAP:])
            p[MWID-RIGHT_OLAP-1,:] = (p[MWID-RIGHT_OLAP-1,:] + np.array(pwms[i+1][0]))/float(2)
            pwm += p.tolist()            

    return np.array(pwm)

def main():
    mainOutDir = '../my_results/zf-C2H2_250_50_seedFFSdiverse6/'

    # Read data
    trainPWMs, trainCores, edges, edges_hmmPos, aaPosList = \
        getPrecomputedInputs_zfC2H2(rescalePWMs=False,ffsOnly=False,includeB1H=False)

    print edges
    print edges_hmmPos
    print aaPosList

    # Remove examples where PWMs that are too short for the number of domains
    nDoms = {}
    for p in trainCores.keys():
        nDoms[p] = len(trainCores[p])/len(aaPosList)
        if len(trainPWMs[p]) < (MWID-RIGHT_OLAP)*nDoms[p]+RIGHT_OLAP or nDoms[p] < 2:
            del nDoms[p]
            del trainCores[p]
            del trainPWMs[p]
    
    # Remove examples where the known/stated fixed starting position
    # would make the pwm too short for the number of arrays annotated as binding
    knownStarts_ffs = readSeedAlignment(SEED_FILE, include = trainPWMs.keys())
    for p in knownStarts_ffs.keys():
        if len(trainPWMs[p]) < knownStarts_ffs[p]['start']+(MWID-RIGHT_OLAP)*nDoms[p]+RIGHT_OLAP:
            #print p, core[p], nDoms[p], len(pwms[p])
            del nDoms[p]
            del trainCores[p]
            del trainPWMs[p]
            del knownStarts_ffs[p]
    fixedStarts = knownStarts_ffs
    tmp = fixedStarts
    fixedStarts = {k: tmp[k] for k in DIVERSE_ZF_SET}

    # Assign to observation groups with identical core sequences
    obsGrps = assignObsGrps(trainCores, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots

    print fixedStarts
    #"""
    # Generate scores for models with random starting points
    fout = open(mainOutDir+'randomModelScores.txt', 'w')
    fout.write('score\n')
    fullX, grpInd = formGLM_fullX(trainCores, edges, uniqueProteins, obsGrps)
    for k in range(250):
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
    #"""

if __name__ == '__main__':
    main()