# A script for making predictions based on parameter
# estimates output by gibbAlign_GLM.py using proteins from CIS-BP

import numpy as np
import os, sys
from pwm import makeNucMatFile, makeLogo
from gibbsAlign_GLM import getHomeoboxData, makeAllLogos
from gibbsAlign_GLM import getAlignedPWMs, getOrientedPWMs, reverseBaseOrient
from getHomeoboxConstructs import subsetDict
from scipy.stats import pearsonr
from copy import deepcopy
from gibbsAlign_GLM import formGLM_fullX, formGLM_testX
from matAlignLib import comp_matrices, matrix_compl, PCC
from scipy.stats import entropy
import pickle
import scipy.stats
import sklearn.metrics

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
TST_SET = 'b08'  # Just controls where code looks for the PCMs

MWID = 6  # Number of base positions in the PDSIM
RAND_SEED = 382738375 #78374223 # Random seed for reproducibilty
N_CHAINS = 200
OUTPUT_LOGO = False


# Helper functions
def getClosestCisProt(b08core, ciscore, cispwms):
    """ Returns a dictionary mapping labels from dict seq to
    a (label, distance) pair corresponding to an entry from
    the ref dict that is closest in terms of hamming distance
    """
    closest = {}
    b08key_cispwms = {}
    b08key_ciscore = {}
    for k in b08core.keys():
        minDist = 1e9
        minKey = ''
        x = b08core[k]
        for r in ciscore.keys():
            dist = 0
            y = ciscore[r]
            for i in range(len(x)):
                if x[i] != y[i]:
                    dist += 1
            if dist < minDist:
                minKey = r
                minDist = dist
        closest[k] = (minKey, minDist)
        b08key_cispwms[k] = cispwms[minKey]
        b08key_ciscore[k] = ciscore[minKey]

    return closest, b08key_cispwms, b08key_ciscore


def bestAliScoreWithNoyes(closestCore, aliPWMS, n08pwms, logoDir, createLogo=False):
    # What is the average per-column pearson correlation
    # for the best IC-PCC alignment to the fly pwm?
    pccVal = {}
    if createLogo == True:
        print("Creating motif logos for closest mouse and cis-bp proteins ...")
    for idx, k1 in enumerate(sorted(closestCore.keys())):#[:1]:
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

        if createLogo == True:
            logoLab = '_'.join([str(idx), k1])
            makeNucMatFile('./tmp/','tmp',aliPWMS[k1])
            makeLogo('./tmp/tmp.txt',logoDir+logoLab+'_mouse.pdf',
                    alpha = 'dna', colScheme = 'classic')
            os.system('rm %s' %'./tmp/tmp.txt')

            logoLab2 = '_'.join([str(idx), k2])
            makeNucMatFile('./tmp/','tmp', n08pwms[k2])
            makeLogo('./tmp/tmp.txt',logoDir+logoLab2+'_cisbp.pdf',
                    alpha = 'dna', colScheme = 'classic')
            os.system('rm %s' %'./tmp/tmp.txt')

    return pccVal

def createLogos(cispwms, pred_pwms, logoDir):
    for p in cispwms.keys():
        makeNucMatFile('./tmp/','tmp',cispwms[p])
        makeLogo('./tmp/tmp.txt',logoDir+p+'_cisbp.pdf',
                alpha = 'dna', colScheme = 'classic')
        os.system('rm %s' %'./tmp/tmp.txt')

        makeNucMatFile('./tmp/','tmp', pred_pwms[p])
        makeLogo('./tmp/tmp.txt',logoDir+p+'_predicted.pdf',
                alpha = 'dna', colScheme = 'classic')
        os.system('rm %s' %'./tmp/tmp.txt')



def compute_alignmentscore(pred_pwm, cisbp_pwm, p):
    s_len = cisbp_pwm.shape[0] - pred_pwm.shape[0] + 1
    score = {}
    idx = 0
    for o in range(2):
        if o == 1:
            cisbp_pwm = matrix_compl(cisbp_pwm)
        for s in range(s_len):
            score[idx] = 0
            for i in range(MWID):
                score[idx] += 0.5 * (2 - entropy(cisbp_pwm[i+s,],base=2)) * PCC(pred_pwm[i,], cisbp_pwm[i+s,])
            idx += 1
    maxs = max(score, key=score.get)
    maxo = 0
    if maxs >= s_len:
        maxo = 1
    if s_len == 0:
        maxs = 0
    else:
        maxs = maxs%s_len
    return maxs, maxo



def main():

    # Obtain Model
    filename = '../results/glmres_chain200.pickle'
    with open(filename) as f:
        res = pickle.load(f)
    score = [x['score'] for x in res]
    opt = np.argmax(score)
    reorient = [x['reorient'] for x in res][opt]
    start = [x['start'] for x in res][opt]
    rev = [x['rev'] for x in res][opt]
    model = [x['final_model'] for x in res][opt]


    # Read data
    if CORE_POS == 'canon9':
        aaPosList = CANON9
    else:
        aaPosList = CORE_POS
    seqs, pwms, core, full, trunc, aaPosList, edges, edges_hmmPos = \
        getHomeoboxData(['n08', 'b08', 'cisbp'], MWID, aaPosList = aaPosList,
                        aposCut = APOS_CUT, edgeCut = EDGE_CUT,
                        maxEdgesPerBase = MAX_EDGES_PER_BASE,
                        N51A_bpos = (MWID-6)/2+2)


    # Extract data for cisbp, subset
    cisseqs, cispwms, ciscore = seqs['cisbp'], pwms['cisbp'], core['cisbp']
    b08seqs, b08pwms, b08core = seqs['b08'], pwms['b08'], core['b08']

    # Discard proteins with amino acid `X` in the core sequence for cisbp and b08
    keep = set([k for k in ciscore.keys() if 'X' not in ciscore[k]])
    [subsetDict(x, keep) for x in [cispwms, ciscore]]
    # Discard proteins with amino acid `-` in the core sequence for cisbp
    keep1 = set([k for k in ciscore.keys() if '-' not in ciscore[k]])
    [subsetDict(x, keep1) for x in [cispwms, ciscore]]
    keep2 = set([k for k in cispwms.keys() if cispwms[k].shape[0] > MWID])
    [subsetDict(x, keep2) for x in [cispwms, ciscore]]
    # Discard cis-bp proteins that are already in mouse proteins (b08)
    keep3 = set([k for k in ciscore.keys() if k not in b08core.keys()])
    [subsetDict(x, keep3) for x in [cispwms, ciscore]]

    # Discard proteins with amino acid `X` in the core sequence for cisbp and b08
    retain = set([k for k in b08core.keys() if 'X' not in b08core[k]])
    [subsetDict(x, retain) for x in [b08pwms, b08core]]
    # Discard proteins with amino acid `-` in the core sequence for cisbp
    retain1 = set([k for k in b08core.keys() if '-' not in b08core[k]])
    [subsetDict(x, retain1) for x in [b08pwms, b08core]]
    retain2 = set([k for k in b08pwms.keys() if b08pwms[k].shape[0] > MWID])
    [subsetDict(x, retain2) for x in [b08pwms, b08core]]

    fullX = formGLM_fullX(ciscore, edges, cispwms.keys())
    pred_pwms = {}
    for idx, p in enumerate(cispwms.keys()):
        pwm = []
        testX = formGLM_testX(fullX, idx)
        for j in range(MWID):
            prediction = model[j].predict_proba(testX[j])[0].tolist()
            pwm.append(prediction)
        pred_pwms[p] = np.array(pwm)

    start, rev = {}, {}
    for p in cispwms.keys():
        start[p], rev[p] = compute_alignmentscore(pred_pwms[p], cispwms[p], p)

    alignedcispwms = getAlignedPWMs(getOrientedPWMs(cispwms, rev), ciscore, start, MWID, flipAli = False)

    pcc_agree, pcc_data, mse_data = {}, {}, {}
    wrongprot_aa = {}
    trainprot_aa = {}
    for j in range(MWID):
        match = 0
        pcc_data[j], mse_data[j] = [], []
        wrongprot_aa[j] = {}
        trainprot_aa[j] = {}
        for aapos in edges[j]:
            wrongprot_aa[j][aapos] = []
            trainprot_aa[j][aapos] = []

        for p in cispwms.keys():
            pcc_val = scipy.stats.pearsonr(pred_pwms[p][j,], alignedcispwms[p][j,])[0]
            pcc_data[j].append(pcc_val)
            mse_val = sklearn.metrics.mean_squared_error(alignedcispwms[p][j,], pred_pwms[p][j,])
            mse_data[j].append(mse_val)
            if pcc_val >= 0.5:
                match += 1
            else:
                for aapos in edges[j]:
                    wrongprot_aa[j][aapos].append(ciscore[p][aapos])

        for p in b08pwms.keys():
            for aapos in edges[j]:
                trainprot_aa[j][aapos].append(b08core[p][aapos])

        for aapos in edges[j]:
            wrongprot_aa[j][aapos] = set(wrongprot_aa[j][aapos])
            trainprot_aa[j][aapos] = set(trainprot_aa[j][aapos])
        pcc_agree[j] = match*0.1 / len(cispwms.keys()) * 10

    aa_compare = {}
    i = 0
    for dic in wrongprot_aa.values():
        for aa_set in dic.values():
            aa_compare[i] = [A2IND[x] for x in list(aa_set)]
            i += 1

    i = 0
    for dic in trainprot_aa.values():
        for aa_set in dic.values():
            aa_compare[i].append([A2IND[x] for x in list(aa_set)])
            i += 1






    # Print results
    print("We are predicting ", len(cispwms.keys()), " proteins.")
    print("pcc_data", pcc_data)
    print("mse_data", mse_data)
    print("pcc_agree", pcc_agree)
    print("edges: ", edges)
    print("wrong_protein aa", wrongprot_aa)
    print("train_protein aa", trainprot_aa)
    print("aa_compare", aa_compare)

    # Save data for analysis
    analysisDir = "../cis_bp/analysis_stats/"
    if not os.path.exists(analysisDir):
        os.makedirs(analysisDir)
    np.savetxt(analysisDir+"pcc_agree_proportion.out", pcc_agree.values())
    for j in range(MWID):
        np.savetxt(analysisDir+"pcc_data_pos"+ str(j) + ".out", pcc_data[j])
        np.savetxt(analysisDir+"mse_data_pos"+ str(j) + ".out", mse_data[j])


    # Create motif logos
    if OUTPUT_LOGO == True:
        logoDir = '../cis_bp/pairlogos/'
        if os.path.exists(logoDir):
            os.rmdir(logoDir)
        os.makedirs(logoDir)
        createLogos(alignedcispwms, pred_pwms, logoDir)



if __name__ == '__main__':
    main()



















