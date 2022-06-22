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
from gibbsAlign_GLM import formGLM_trainX, formGLM_trainW, formGLM_Y, createGLMModel
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

def predictSpecificity_array_ZF(fullX, model, startInd, arrayLen, wtB1 = 0.5):
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

    # Weighted average for overlapping positions for adjacent domains and concatenate
    pwm = []
    for i in range(arrayLen):
        if i == 0:
            p = np.array(pwms[i])
            p[MWID-1,:] = (1-wtB1)*p[MWID-RIGHT_OLAP,:] + wtB1*np.array(pwms[i+1][0])
            pwm += p.tolist()            
        elif i == arrayLen-1:
            pwm += pwms[i][RIGHT_OLAP:]
        else:
            p = np.array(pwms[i][RIGHT_OLAP:])
            p[MWID-RIGHT_OLAP-1,:] = (1-wtB1)*p[MWID-RIGHT_OLAP-1,:] + wtB1*np.array(pwms[i+1][0])
            pwm += p.tolist()            

    return np.array(pwm)

def main():
    mainOutDir = '../my_results/zf-C2H2_250_50_seedFFSdiverse6_80_20/'

    # Obtain Model
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
    trainPWMs, trainCores, edges, edges_hmmPos, aaPosList = \
        getPrecomputedInputs_zfC2H2(rescalePWMs=False,ffsOnly=False,includeB1H=False)

    print edges
    print edges_hmmPos
    print aaPosList
    #print trainCores

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

    # Assign to observation groups with identical core sequences
    obsGrps = assignObsGrps(trainCores, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots

    # Create the predicted logos for the corresponding core sequences
    # while holding out the core seq grps
    print("Creating predicted PWMs according to retrained models ...")
    fullX, grpInd = formGLM_fullX(trainCores, edges, uniqueProteins, obsGrps, domainOrder = -1)

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

        # If anchoring on the B1H data, skip that data in HOO vaildation testing
        if len(testProteins) == 1 and ANCHOR_B1H and testProteins[0][:4] == 'B1H.':
                continue

        trainX = formGLM_trainX(fullX,startInd_ho,endIndex_ho)
        trainY = formGLM_Y(trainProteins, nDoms)
        trainW = formGLM_trainW(trainPWMs, trainProteins, nDoms, start, rev)
        model_ho = createGLMModel(trainX, trainY, trainW)
        
        #print coreSeq
        pwm = predictSpecificity_array_ZF(fullX, model_ho, startInd_ho, 
                                          nDoms[testProteins[0]], wtB1 = 0.5)

        for p in testProteins:
            pred_pwms[p] = pwm
            pWid = len(pwm)
            # Align to the predicted PWM
            score, shift, ori = comp_matrices(pred_pwms[p], trainPWMs[p],
                                              minWidthM2 = pWid, oneSided=True)
            shift = -shift
            if ori:
                test_pwms_ali[p] = matrix_compl(trainPWMs[p])[shift:shift+pWid]
            else:
                test_pwms_ali[p] = trainPWMs[p][shift:shift+pWid]
            assert len(pred_pwms[p]) == len(test_pwms_ali[p])

    print("Ran in %.2f seconds" %(time.time()-startTime))
    pccTabfile = mainOutDir+'/pccTable_underS_holdOneOut.txt'
    makePCCtable(test_pwms_ali, pred_pwms, trainCores, pccTabfile)
    makePWMtab(pred_pwms,mainOutDir+'/pwms_pred_holdOneOut.txt')
    makePWMtab(test_pwms_ali,mainOutDir+'/pwms_testAli_holdOneOut.txt')

    print len(test_pwms_ali), len(pred_pwms)
    # Make predicted logos from the hold-one-out training
    if MAKE_LOGOS:
        print "Generating the logos from hold-one-out training"
        logoDir = mainOutDir + '/0_logos_predicted_hoo/'
        logoDir2 = mainOutDir + '/0_logos_aligned_hoo/'
        startTime = time.time()
        if not os.path.exists(logoDir):
            os.makedirs(logoDir)
        makeAllLogos(pred_pwms, trainCores, logoDir)
        makeAllLogos(test_pwms_ali, trainCores, logoDir2)
        print("Ran in %.2f seconds" %(time.time()-startTime))

if __name__ == '__main__':
    main()

