# Apply specififcity transfer for PWMs we have alignments for 
# against PWMs of proteins that are HD1 away in terms of core sequence

from matAlignLib import *
import numpy as np
import scipy
import os, sys, re
from pwm import makeNucMatFile, makeLogo
from gibbsAlign_GLM import getHomeoboxData, makeCoefTable, makePWMtab
from gibbsAlign_GLM import getAlignedPWMs, getOrientedPWMs
from gibbsAlign_GLM import assignObsGrps, computeGLMLoglikelihood
from gibbsAlign_GLM import formGLM_fullX, formGLM_testX, formGLM_testW
from gibbsAlign_GLM import getTrainPairsAndInfo, getTestProtsAndPWMs
from getHomeoboxConstructs import subsetDict
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
COMBINE = 'fracID_fullDBD'#'mean'  #
MAKE_LOGOS = True
RESCALE_PWMS = True

def readHD1proteinTable(fname):
    # Returns a dictionary contiaining information about test proteins
    # and ther HD-1 neighbors in the training set

    fin = open(fname, 'r')
    fin.readline()
    info = {}
    for line in fin:
        testProt, trainProt, subpos, aaFrom, aaTo = line.strip().split('\t')
        subpos = int(subpos.split('.')[-1])
        if info.has_key(testProt):
            if info[testProt].has_key(subpos):
                info[testProt][subpos].append([trainProt, aaTo])
            else:
                info[testProt][subpos] = [[trainProt, aaTo]]
        else:
            info[testProt] = {}
            info[testProt][subpos] = [[trainProt, aaTo]]
    fin.close()
    return info

def hammingDist(a, b):
    hd = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            hd += 1
    return hd

def getNbrWts(toSub, combine = 'mean', dbd_tst_prot = None, 
              dbd_tr = None, maxNbrs = 1e10, minNbrWt = 0.8):
    # Returns a set of relative weights for combining neighbors

    nbrWts = {}
    maxWt = 0
    minWt = 100
    for prot in toSub:

        if combine == 'mean':
            nbrWts[prot] = 1.0
        elif combine == 'fracID_fullDBD':
            #print hammingDist(dbd_tst_prot,dbd_tr[prot])
            nwt = float(len(dbd_tst_prot) - hammingDist(dbd_tst_prot,dbd_tr[prot]))/len(dbd_tst_prot)
            if nwt >= minNbrWt:
                nbrWts[prot] = nwt
                if nwt > maxWt:
                    maxWt = nwt
                if nwt < minWt:
                    minWt = nwt
    #if maxWt < 0.8:
    #    nbrWts = {}
    #for k in nbrWts.keys():
    #    if (maxWt == minWt):
    #        nbrWts[k] = 1.0
    #    else:
    #        nbrWts[k] = (nbrWts[k]-minWt)/(maxWt-minWt)
    return nbrWts

def getNoContactBasePositions(edges, aaPosList):
    # Get mapping from canonical aa positions to base positions NOT CONTACTED
    # by that aa position according to the model's edgeList

    bPosMap = {}
    for bpos in edges.keys():
        for apos in edges[bpos]:
            if bPosMap.has_key(aaPosList[apos]):
                bPosMap[aaPosList[apos]].append(bpos)
            else:
                bPosMap[aaPosList[apos]] = [bpos]
    for aapos in bPosMap.keys():
        bPosMap[aapos] = set(range(MWID)) - set(bPosMap[aapos])
    return bPosMap   

def makeHD1Substitutions(hd1Info, aliPWMs, predPWMs, aaPosList, edges, 
                         combine = 'mean', dbd_tst = None, dbd_tr = None,
                         maxNbrs = 1e10, minNbrWt = 0.80):
    # For each test protein that is HD-1 away from a training protein with 
    # an aligned PWM, save predicted columns only for the base positions 
    # contacted by the HD-1 substituted amino acid position, then replace the
    # remaining base positions with a combination of the PWM columns for 
    # that aligned position across all of the HD-1 neighbors of the test protein.
    # Also returns the list of proteins used when subsituting for each base position.

    bPosMap = getNoContactBasePositions(edges, aaPosList)

    # Make the transfers for each predPWMs according to the bPosMap
    # I.e., for each of base position for a test protein, prefer to 
    # make a substituion based on the set of aligned HD-1 neighbor PWMs
    # for proteins whose aa substitution positions do not contact
    # that base positon according to the model's edgeList.
    subPWMs = {}
    nbrWts_all = {}
    for prot in hd1Info.keys():#[:1]:
        nbrWts_all[prot] = {}
        toSub = {}
        for bpos in range(MWID):
            toSub[bpos] = []
            for aapos in hd1Info[prot].keys():
                if bpos in bPosMap[aapos]:
                    toSub[bpos] += [x[0] for x in hd1Info[prot][aapos]]
        #for bpos in sorted(toSub.keys()):
        #    print prot, bpos, toSub[bpos]

        #if prot == '11_TGA_KSVAQ':
        #    print toSub

        new_pwm = np.zeros((MWID,4), dtype = 'float')
        for i in range(MWID):
            #print i
            if len(toSub[bpos]) == 0:
                new_pwm[i] = predPWMs[prot][i]
            else:
                nbrWts = getNbrWts(toSub[i], combine = combine, 
                                   dbd_tst_prot = dbd_tst[prot], 
                                   dbd_tr = dbd_tr, maxNbrs = len(toSub[bpos]),
                                   minNbrWt = minNbrWt)
                #print prot, i, len(nbrWts.keys())
                nbrWts_all[prot][i] = nbrWts
                #print prot, i, len(nbrWts_all[prot][i].keys())
                #print nbrWts
                for nbr in nbrWts:
                    if nbrWts.has_key(nbr):
                        new_pwm[i] += aliPWMs[nbr][i]*nbrWts[nbr]
            
            # Correct for case where no neighbors had enough identity
            if np.all(new_pwm[i] == 0):
                new_pwm[i] = predPWMs[prot][i]

            new_pwm[i] = new_pwm[i]/new_pwm[i].sum()
        subPWMs[prot] = new_pwm
    return subPWMs, nbrWts_all

def makeAllLogos(pwms, logoDir, keysToUse = None):
    """
    Output function
    Place logos for every pwm and ouput to logoDir
    """

    if keysToUse == None:
        keysToUse = sorted(pwms.keys())

    if not os.path.exists(logoDir):
        os.makedirs(logoDir)

    i = 0
    for k in keysToUse:
        pwm = pwms[k]
        #if k == 'PBX1':
        #    print pwm
        logoLab = k+'_%d' %i
        makeNucMatFile('./tmp/','tmp',pwm)
        makeLogo('./tmp/tmp.txt',logoDir+logoLab+'.pdf',
                 alpha = 'dna', colScheme = 'classic')
        os.system('rm %s' %'./tmp/tmp.txt')
        i += 1

def makePCCtable(exp, pred, fname, nbrWts = None, 
                 transCols = None, aaPosList = None):
    # Compute the per-position pearson correlation coefficients
    # and mean squared errors between the predicted and the 
    # experimental pwms assuming that the alignment is correct

    fout = open(fname, 'w')
    fout.write('prot\tpos\tpcc\trmse\tic.exp\tic.pred\tnnbrs')
    if transCols is not None:
        fout.write('\ttransfer\tvarPos\n')
    else:
        fout.write('\n')
    for prot in sorted(exp.keys()):
        for pos in range(len(exp[prot])):
            x, y = exp[prot][pos], pred[prot][pos]
            pcc = scipy.stats.pearsonr(x,y)[0]
            rmse = np.sqrt(np.mean((x - y)**2))
            icExp, icPred = information_content(x), information_content(y)
            fout.write('%s\t%d\t%f\t%f\t%f\t%f' %(prot, pos+1, pcc, rmse, icExp, icPred))

            if nbrWts is not None and nbrWts.has_key(prot) and nbrWts[prot].has_key(pos):
                nnbrs = len(nbrWts[prot][pos].keys())
            else:
                nnbrs = 0

            if transCols is None:
                fout.write('\t%d\n' %(nnbrs))
            else:
                if pos in transCols[prot]['wtCols']:
                    transfer = True
                    varPos = aaPosList[transCols[prot]['varPos']]
                else:
                    transfer = False
                    varPos = aaPosList[transCols[prot]['varPos']]
                fout.write('\t%d\t%s\t%d\n' %(nnbrs,transfer,varPos))
    fout.close()

def readStormoPWMs(fnames):
    # Returns a set of numpy matrices for PWMs from a list of 
    # text files each with output in the format of the RF models 
    # from Gary Stormo's group

    pwms = {}
    for fname in fnames:
        fin = open(fname, 'r')
        line = fin.readline()
        while True:
            if line == '':
                break
            elif line[0] == '#':
                line = fin.readline()
                continue
            elif line[0] == '>':
                pname = line[1:].split('/')[0]
                a = [float(x) for x in fin.readline().split('|')[-1].strip().split()]
                c = [float(x) for x in fin.readline().split('|')[-1].strip().split()]
                g = [float(x) for x in fin.readline().split('|')[-1].strip().split()]
                t = [float(x) for x in fin.readline().split('|')[-1].strip().split()]
                pwms[pname] = np.zeros((len(a),4), dtype = float)
                for i in range(len(a)):
                    pwms[pname][i] = [a[i],c[i],g[i],t[i]]
                pwms[pname]
                line = fin.readline()
            else:
                line = fin.readline()
                continue
    return pwms

def compareTestToStormoPreds(pwms_test, predDir):
    # Aligns and makes comparisons of predictions from Gary Stormo's
    # lab to the test set PWMs

    fnames = [predDir + x for x in os.listdir(predDir)]
    pred_pwms = readStormoPWMs(fnames)
    testLabs = set(pwms_test.keys())
 
    pwms_test_pred = {}    
    # Align to the test PWM and make comparisons
    for p in pred_pwms.keys():
        
        # Since some could not be predicted by RF
        if p not in testLabs:
            #pwms_test_pred[p] = None
            continue

        pred_pwms[p] = pred_pwms[p][2:2+MWID]
        score, shift, ori = comp_matrices(pred_pwms[p], pwms_test[p],
                                          minWidthM2 = 6, oneSided=True)
        shift = -shift
        if ori:
            pwms_test_pred[p] = matrix_compl(pwms_test[p])[shift:shift+MWID]
        else:
            pwms_test_pred[p] = pwms_test[p][shift:shift+MWID]   

    makePCCtable(pwms_test_pred, pred_pwms, predDir+'pccTable_test.txt')

def findMaxStartPos(prot, pwm, cs, edges, model):
    # Find the maximum starting position for a pwm given the model and 
    # it core sequence

    seqWid = len(pwm)
    mWid = len(edges.keys())
    pwmDict = {prot: pwm}
    coreDict = {prot: cs}
    obsGrps = assignObsGrps(coreDict)
    fullX, grpInd = formGLM_fullX(coreDict, edges, [prot], obsGrps)

    maxll, maxs, maxr = -1e30, 0, 0
    #print pwm
    for r in range(2):
        for s in range(seqWid-mWid+1):
            start = {prot: s}
            rev = {prot :r}
            testX = formGLM_testX(fullX, 0)
            testW = formGLM_testW(pwmDict, 0, [prot], start, rev)
            ll = computeGLMLoglikelihood(testX, testW, model)
            #print r, s, ll, maxll
            if ll > maxll:
                maxll = ll
                maxs = s
                maxr = r
    return maxs, maxr

def tranferPredBarerra2012(pwms_test, pwms_train_aligned, core_test, 
                           core_train, edges, model, subPWMs, pred_pwms):
    # What if, for each Barerra mutant (HD-1, varying in a soecificty residue), 
    # we just take each column from the wildtype version of that protein
    # except for those contacting then do either nn transfer or model prediction
    # for the remaning columns' predictions.
    # 
    # If the wildtype protein was held out of the training set due to 
    # having identical core sequence to another protein in the test set
    # (but not the mutant), then we take the alignment that mazimizes
    # its likelihood to be the best one (i.e., same procedure as sampling
    # to align, but just set the max prob position to be 1)

    # Find the barerra mutant proteins using regex
    p = '.*_[A-Z][0-9]{1,3}[A-Z]$'
    #print pwms_test.keys()
    #print re.match(p, 'HESX1_N125S')
    muts = [x for x in pwms_test.keys() if re.match(p,x) is not None]
    print(len(muts))
    wts = [x.split('_')[0] for x in muts]

    # Find the varying core position, if there is one, as well as 
    # whether the wild-type protein is in the training or testing set
    varPos = {}
    inTest = {}
    for i in range(len(muts)):
        wt, mut = wts[i], muts[i]
        if wt in core_test.keys():
            a = core_test[wt]
            inTest[wt] = True
        elif wt in core_train.keys():
            a = core_train[wt]
            inTest[wt] = False
        else:
            pass
        b = core_test[mut]
        for i in range(len(a)):
            if a[i] != b[i]:
                assert (wt,mut) not in varPos.keys()
                varPos[(wt,mut)] = i

    # Find the optimal alignments for the wt proteins
    aliWT = {}
    for k in varPos.keys():
        wt = k[0]
        if inTest[wt]:
            cs = core_test[wt]
            s, r = findMaxStartPos(wt, pwms_test[wt], cs, edges, model)
            if r:
                wtPWM = matrix_compl(pwms_test[wt])[s:s+MWID]
            else:
                wtPWM = pwms_test[wt][s:s+MWID]
        else:
            wtPWM = pwms_train_aligned[wt]
        aliWT[wt] = wtPWM

    # For base positions that are not contacted by the substituted amino acid 
    # position, transfer the corresponding WT column to that position
    # in each of the predicted mutant PWMs
    # (for both model-based prediction and for model+nn)
    aaPosList = range(len(a))
    bPosMap = getNoContactBasePositions(edges, aaPosList)
    transfer_model = {}
    transfer_modelNN = {}
    wtCols = {}
    for (wt, mut) in varPos.keys():
        if subPWMs.has_key(mut):
            transfer_modelNN[mut] = np.copy(subPWMs[mut])
        else:
            transfer_modelNN[mut] = np.copy(pred_pwms[mut])
        transfer_model[mut] = np.copy(pred_pwms[mut])
        wtCols[mut] = {'wtCols': list(bPosMap[varPos[(wt, mut)]]),
                       'varPos': varPos[(wt,mut)]}
        for i in list(bPosMap[varPos[(wt, mut)]]):
            transfer_modelNN[mut][i] = aliWT[wt][i]
            transfer_model[mut][i] = aliWT[wt][i]

    return transfer_model, transfer_modelNN, wtCols

def tranferAllHD1Pairs(hd1Info, pwms_test, pwms_train_aligned, core_test,
                       core_train, full_test, full_train, edges, 
                       aaPosList, model, pred_pwms, simThresh = 0.8):
    # For each pair of proteins within a given overall DBD similarity
    # threshold, varying in exactly one specificity residue, and on opposite
    # sides of the train/test partition, transfer base positions that don't 
    # touch the varied specificity residue, then predict for the remaining 
    # base position using the model or the NN routine

    bPosMap = getNoContactBasePositions(edges, aaPosList)   
    transferPWMs = {}
    for prot1 in hd1Info.keys():
        transferPWMs[prot1] = {}
        for varPos in hd1Info[prot1].keys():
            transferPWMs[prot1][varPos] = {}
            replace = set(range(MWID)) - bPosMap[varPos]
            for (prot2, aa) in hd1Info[prot1][varPos]:
                #print prot1, varPos, prot2, aa
                hDist_full = hammingDist(full_test[prot1], full_train[prot2])
                dbdLen = len(full_test[prot1])
                dbdSim = (dbdLen-hDist_full)/float(dbdLen)
                if dbdSim < simThresh:
                    continue
                pwm = np.copy(pwms_train_aligned[prot2])
                for i in replace:
                    pwm[i] = pred_pwms[prot1][i]
                transferPWMs[prot1][varPos][(prot2,aa)] = pwm
    return transferPWMs,bPosMap

def makeTranferAllHD1Pairs_cmpTab(transferPWMs, pwms_test, bPosMap,
                                  pccTabName):
    # Align the transfer PWMs to the test pwms and make per-column 
    # comparisons

    fout = open(pccTabName, 'w')
    fout.write('prot1\tprot2\tpos\tpcc\trmse\tic.exp\tic.pred\ttransfer\tvarPos\taa\n')
    pwms_test_pred = {}
    for prot1 in transferPWMs.keys():
        pwms_test_pred[prot1] = {}
        for varPos in transferPWMs[prot1].keys():
            pwms_test_pred[prot1][varPos] = {}
            for (prot2, aa) in transferPWMs[prot1][varPos].keys():
                # Align to the predicted PWM (model only)
                pwm_pred = transferPWMs[prot1][varPos][(prot2,aa)]
                pwm_test = pwms_test[prot1]
                #print pwm_test
                score, shift, ori = comp_matrices(pwm_pred, pwm_test,
                                                  minWidthM2 = 6, oneSided=True)
                shift = -shift
                if ori:
                    pwmAli = matrix_compl(pwm_test)[shift:shift+MWID]
                else:
                    pwmAli = pwm_test[shift:shift+MWID]
                for i in range(len(pwmAli)):
                    x, y = pwmAli[i], pwm_pred[i]
                    pcc = scipy.stats.pearsonr(x,y)[0]
                    rmse = np.sqrt(np.mean((x - y)**2))
                    icExp, icPred = information_content(x), information_content(y)
                    xfer = False
                    if i in bPosMap[varPos]:
                        xfer = True
                    fout.write('%s\t%s\t%d\t%f\t%f\t%f\t%f\t%s\t%d\t%s\n' \
                               %(prot1, prot2, i+1, pcc, rmse, icExp, 
                                 icPred, xfer, varPos, aa))
                    pwms_test_pred[prot1][varPos][(prot2, aa)] = pwmAli
    return pwms_test_pred


def main():

    mainOutDir = '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50'

    # Obtain Model
    outLabel = 'cisbp-chu/structFixed1_grpHoldout_multinom_chain100maxIter15scaled50'
    filename = '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15scaled50.pickle'
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

    # Get the oritented and aligned 6bp pwms from the training set
    # and information about the edges in the model and positions 
    # in the core sequences
    pwms_train, core_train, full_train, edges, edges_hmmPos, aaPosList, testProts = \
        getTrainPairsAndInfo(rescalePWMs = RESCALE_PWMS)
    print len(pwms_train), len(start)
    assert len(pwms_train) == len(start)
    pwms_train_oriented = getOrientedPWMs(pwms_train, rev)
    pwms_train_aligned = getAlignedPWMs(pwms_train_oriented,core_train, 
                                        start, MWID, flipAli = False)

    # Get the ground truth test proteins and PWMs
    pwms_test, core_test, full_test = getTestProtsAndPWMs(testProts)

    # Make predicted PWMs based on each of the test core sequences
    obsGrps = assignObsGrps(core_test, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots 
    fullX, grpInd = formGLM_fullX(core_test, edges, uniqueProteins, obsGrps)
    pred_pwms = {}
    for idx, p in enumerate(uniqueProteins):
        pwm = []
        testX = formGLM_testX(fullX, idx)
        for j in range(MWID):
            prediction = model[j].predict_proba(testX[j])[0].tolist()
            pwm.append(prediction)
        pred_pwms[p] = np.array(pwm)

    # Read in the names of the HD-1 test protein information
    hd1Info = readHD1proteinTable('../hd1-explore/0_splitChu2012-trainTest/0_hd1cores_test_train.txt')
    assert len(set(hd1Info.keys()) - set(pwms_test.keys())) == 0

    # Make the appropriate substitutions into the predicted PWMs
    # based on the HD-1 neighbors of the test protein in the training set
    subPWMs, nbrWts = \
        makeHD1Substitutions(hd1Info, pwms_train_aligned, pred_pwms, 
                             aaPosList, edges, combine = COMBINE,
                             dbd_tst = full_test, dbd_tr = full_train,
                             minNbrWt = 0.80)
    
    # Find the optimal alignment between the substitution PWM and the 
    # corresponding ground truth PWM and measure accuracy
    pwms_test_sub = {}   # test pwms aligned to substituiton pwms
    pwms_test_pred = {}  # test pwms aligned to predicted pwms
    for p in subPWMs.keys():

        # Align to the transfer PWM
        score, shift, ori = comp_matrices(subPWMs[p], pwms_test[p],
                                          minWidthM2 = 6, oneSided=True)
        shift = -shift
        if ori:
            pwms_test_sub[p] = matrix_compl(pwms_test[p])[shift:shift+MWID]
        else:
            pwms_test_sub[p] = pwms_test[p][shift:shift+MWID]

        # Align to the predicted PWM (model only)
        score, shift, ori = comp_matrices(pred_pwms[p], pwms_test[p],
                                          minWidthM2 = 6, oneSided=True)
        shift = -shift
        if ori:
            pwms_test_pred[p] = matrix_compl(pwms_test[p])[shift:shift+MWID]
        else:
            pwms_test_pred[p] = pwms_test[p][shift:shift+MWID]
        #if p == '11_TGA_KSVAQ':
        #    print shift, ori
        assert len(pwms_test_sub[p]) == MWID
        assert len(pwms_test_pred[p]) == MWID

    outDir = mainOutDir+'/transfer_test/%s/'%COMBINE
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    makePCCtable(pwms_test_sub, subPWMs, outDir+'pccTable_test_transfer.txt', 
                 nbrWts = nbrWts)
    makePCCtable(pwms_test_pred, pred_pwms, outDir+'pccTable_test_predOnly.txt',
                 nbrWts = nbrWts)
    makePWMtab(subPWMs,outDir+'/pwms_glmNN_testSet.txt')
    makePWMtab(pred_pwms,outDir+'/pwms_glmOnly_testSet.txt')

    # Output a table for exploring the model coefficients
    #print model
    makeCoefTable(model, edges_hmmPos, outDir+'modelCoefs.txt')

    #print pwms_test_sub['PBX1']
    #print pwms_test_pred['PBX1']
    
    if MAKE_LOGOS:
        makeAllLogos(pwms_test_sub,outDir+'transfer_testLogos_bestAli/')
        makeAllLogos(subPWMs,outDir+'transfer_predLogos_bestAli/')
        makeAllLogos(pwms_test_pred,outDir+'predOnly_testLogos_bestAli/')
        makeAllLogos(pred_pwms,outDir+'predOnly_predLogos_bestAli/')
    #"""

    # Align and compare the predicitons from Gary Stormo's group
    compareTestToStormoPreds(pwms_test,'../stormo_predictor_outputs/extant/')
    compareTestToStormoPreds(pwms_test,'../stormo_predictor_outputs/joint/')

    # See how we do predicting HD-1 with mutant vs. single wildtype prot
    #print subPWMs.keys()
    transfer_model, transfer_modelNN, wtCols = \
        tranferPredBarerra2012(pwms_test, pwms_train_aligned, core_test, 
                               core_train, edges, model, subPWMs, pred_pwms)
    barreraMutPWMs_test_model = {}
    barreraMutPWMs_test_modelNN = {}
    for p in transfer_model.keys():
        
        # Align to the predicted PWM (model only)
        score, shift, ori = comp_matrices(transfer_model[p], pwms_test[p],
                                          minWidthM2 = 6, oneSided=True)
        shift = -shift
        if ori:
            barreraMutPWMs_test_model[p] = matrix_compl(pwms_test[p])[shift:shift+MWID]
        else:
            barreraMutPWMs_test_model[p] = pwms_test[p][shift:shift+MWID]
        
        # Align to the predicted PWM (model + nn)
        score, shift, ori = comp_matrices(transfer_modelNN[p], pwms_test[p],
                                          minWidthM2 = 6, oneSided=True)
        shift = -shift
        if ori:
            barreraMutPWMs_test_modelNN[p] = matrix_compl(pwms_test[p])[shift:shift+MWID]
        else:
            barreraMutPWMs_test_modelNN[p] = pwms_test[p][shift:shift+MWID]    
    print outDir
    makePCCtable(barreraMutPWMs_test_model, transfer_model, 
                 outDir+'pccTable_BarreraMuts_transfer-model.txt',
                 transCols = wtCols, aaPosList = aaPosList)
    makePCCtable(barreraMutPWMs_test_modelNN, transfer_modelNN, 
                 outDir+'pccTable_BarreraMuts_transfer-modelNN.txt',
                 nbrWts = nbrWts,transCols = wtCols, aaPosList = aaPosList)
    #for k in sorted(wtCols.keys()):
    #    print k, wtCols[k]

    # Test transfer mehtod for all HD1 pairs across the train/test partition
    transferHD1, bPosMap = \
        tranferAllHD1Pairs(hd1Info, pwms_test, pwms_train_aligned, core_test,
                           core_train, full_test, full_train, edges, 
                           aaPosList, model, pred_pwms, simThresh = 0.8)
    transferHD1_NN, bPosMap = \
        tranferAllHD1Pairs(hd1Info, pwms_test, pwms_train_aligned, core_test,
                           core_train, full_test, full_train, edges, 
                           aaPosList, model, subPWMs, simThresh = 0.8)
    transferHD1_ali = \
        makeTranferAllHD1Pairs_cmpTab(transferHD1, pwms_test, bPosMap,
                                      outDir+'pccTable_allHD1Pairs_transfer-model.txt')
    transferHD1_NN_ali = \
        makeTranferAllHD1Pairs_cmpTab(transferHD1_NN, pwms_test, bPosMap,
                                      outDir+'pccTable_allHD1Pairs_transfer-modelNN.txt')



if __name__ == '__main__':
    main()

    
