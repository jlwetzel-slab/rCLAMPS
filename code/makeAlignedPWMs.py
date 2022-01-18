# Outputs the aligned experimental PWMs according to the alignment produced 
# by gibbsAlign_GLM.py .  Also outputs corresponding predicted PWMs from M_theta
# based on only the core amino acid sequence residues of the corresponding 
# protiens.  Then computed information regarding the agreement between PWMs 
# predicted by M_theta to the experimetal PWMs, assuming that the alignment 
# produced by gibbsAlign_GLM.py was correct.

# For validating the alignment itself, we use the GLM_predict.py script 
# to compare to the experimentally aligned orthologous fly PWMs from Noyes 08 (Cell).

import numpy as np
import scipy
import os, sys
from pwm import makeNucMatFile, makeLogo
from gibbsAlign_GLM import getHomeoboxData, makeAllLogos, makePWMtab
from gibbsAlign_GLM import getAlignedPWMs, getOrientedPWMs
from gibbsAlign_GLM import assignObsGrps, formGLM_fullX, formGLM_testX
from gibbsAlign_GLM import getTrainPairsAndInfo
from getHomeoboxConstructs import subsetDict
from holdoutRetrainGLM import makePCCtable
import pickle

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

def main():

    if EXCLUDE_TEST:
        dirStem = '../results/cisbp-chu/'
        labStem = 'cisbp-chu/'
    else:
        dirStem = '../results/cisbp-chuAll/'
        labStem = 'cisbp-chuAll/'

    mainOutDir = dirStem+'structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50'

    # Obtain Model
    outLabel = labStem+'structFixed1_grpHoldout_multinom_chain100maxIter15scaled50'
    filename = dirStem+'structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50.pickle'
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
        getTrainPairsAndInfo(rescalePWMs = RESCALE_PWMS, excludeTestSet = EXCLUDE_TEST)
    

    obsGrps = assignObsGrps(trainCores, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots  

    # Output the alignments under the learned vector S alphabetically
    print("Creating aligned experimental logos according to model ...")
    logoDir = mainOutDir + '/experimental_logos_underS/'
    if not os.path.exists(logoDir):
        os.makedirs(logoDir)
    print("Printing logos to: %s" %logoDir)
    flipAli = False
    if reorient:
        flipAli = True
    aliPWMs = getAlignedPWMs(getOrientedPWMs(trainPWMs, rev),
                             trainCores, start, MWID, flipAli = flipAli)
    if MAKE_LOGOS:
        print("Creating experimental logos under alignment vector S ...")
        makeAllLogos(aliPWMs, trainCores, logoDir)

    # Create the predicted logos for the corresponding core sequences
    print("Creating predicted PWMs according to model ...")
    fullX, grpInd = formGLM_fullX(trainCores, edges, uniqueProteins, obsGrps)
    pred_pwms = {}
    for idx, p in enumerate(uniqueProteins):
        pwm = []
        testX = formGLM_testX(fullX, idx)
        for j in range(MWID):
            prediction = model[j].predict_proba(testX[j])[0].tolist()
            pwm.append(prediction)
        pred_pwms[p] = np.array(pwm)
    makePWMtab(pred_pwms,mainOutDir+'/pwms_pred.txt')
    
    # Output the fit information for individual PWMs under S
    pccTabfile = mainOutDir+'/pccTable_underS.txt'
    makePCCtable(aliPWMs, pred_pwms, trainCores, pccTabfile)

    if MAKE_LOGOS:
        print("Creating predicted logos according to model ...")
        logoDir = mainOutDir + '/predicted_logos/'
        if not os.path.exists(logoDir):
            os.makedirs(logoDir)
        makeAllLogos(pred_pwms, trainCores, logoDir)
    


if __name__ == '__main__':
    main()

