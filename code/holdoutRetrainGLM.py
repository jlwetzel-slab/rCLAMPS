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
from gibbsAlign_GLM import getTrainPairsAndInfo
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

    # Create the predicted logos for the corresponding core sequences
    # while holding out the core seq grps
    print("Creating predicted PWMs according to retrained models ...")
    fullX, grpInd = formGLM_fullX(trainCores, edges, uniqueProteins, obsGrps)
    pred_pwms = {}
    test_pwms_ali = {}
    startTime = time.time()
    for k, coreSeq in enumerate(grpInd.keys()):
        startInd_ho, endIndex_ho = grpInd[coreSeq]
        trainProteins, testProteins = [], []
        for prot in uniqueProteins:
            if prot in obsGrps[coreSeq]:
                testProteins += [prot]
            else:
                trainProteins += [prot]
        trainX = formGLM_trainX(fullX,startInd_ho,endIndex_ho)
        trainY = formGLM_Y(trainProteins, nDoms)
        trainW = formGLM_trainW(trainPWMs, trainProteins, nDoms, start, rev)
        model_ho = createGLMModel(trainX, trainY, trainW)
        testX = formGLM_testX(fullX, startInd_ho, startInd_ho + 4)
        pwm = []
        for j in range(MWID):
            prediction = model_ho[j].predict_proba(testX[j])[0].tolist()
            pwm.append(prediction)
        for p in testProteins:
            pred_pwms[p] = np.array(pwm)
            # Align to the predicted PWM
            score, shift, ori = comp_matrices(pred_pwms[p], trainPWMs[p],
                                              minWidthM2 = 6, oneSided=True)
            shift = -shift
            if ori:
                test_pwms_ali[p] = matrix_compl(trainPWMs[p])[shift:shift+MWID]
            else:
                test_pwms_ali[p] = trainPWMs[p][shift:shift+MWID]
    print("Ran in %.2f seconds" %(time.time()-startTime))
    
    # Output the fit information for individual PWMs under S
    if ADD_RESIDUES:
        suff = '_addResCS2012'
    else:
        suff = ''

    pccTabfile = mainOutDir+'/pccTable_underS_holdOneOut'+suff+'.txt'
    makePCCtable(test_pwms_ali, pred_pwms, trainCores, pccTabfile)
    makePWMtab(pred_pwms,mainOutDir+'/pwms_pred_holdOneOut'+suff+'.txt')
    makePWMtab(test_pwms_ali,mainOutDir+'/pwms_testAli_holdOneOut'+suff+'.txt')

    # Make predicted logos from the hold-one-out training
    if MAKE_LOGOS:
        print "Generating the logos from hold-one-out training"
        logoDir = mainOutDir + '/predicted_logos_holdOneOut'+suff+'/'
        startTime = time.time()
        if not os.path.exists(logoDir):
            os.makedirs(logoDir)
        makeAllLogos(pred_pwms, trainCores, logoDir)
        print("Ran in %.2f seconds" %(time.time()-startTime))

if __name__ == '__main__':
    main()
