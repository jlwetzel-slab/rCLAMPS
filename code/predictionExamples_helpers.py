# Helper functions and global variables for predictionExample.py

from sklearn.linear_model import LogisticRegression
from gibbsAlign_GLM import makePWMtab
from getHomeoboxConstructs import readFromFasta, writeToFasta, subsetDict, makeMatchStateTab
from pwm import makeNucMatFile, makeLogo
from matAlignLib import matrix_compl
from runhmmer import runhmmer3, getdescs
import numpy as np
import pickle, os

PROT_SEQ_FILE = '../precomputedInputs/proteins_homeodomains_hasPWM.fa'  # Input protein sequence fasta file
PWM_INPUT_TABLE = '../precomputedInputs/pwmTab_homeodomains_all.txt'   # A table of PWMs corresponding to prots in PROT_SEQ_FILE
CONTACT_MAP = '../precomputedInputs/homeodomain_contactMap.txt'  # A contact map for the domain family
HMM_FILE = '../pfamHMMs/Homeobox.hmm'    # Location of hmm file
HMM_LEN = 57                             # Number of match states in HMM file
HMM_NAME = 'Homeobox'                    # A name for the HMM
HMM_OFFSET = 2                           # Used to offset HMM states to canonical numbering scheme for Homoedomains
HMMER_HOME = None
EXCLUDE_TEST = False
OBS_GRPS = 'grpIDcore'
MWID = 6

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

def readPWMtab(fname):
    # Reads in a PWM table in the format output by makePWMtab
    fin = open(fname, 'r')
    fin.readline()
    pwms = {}
    for line in fin:
        prot, bpos, base, freq = line.strip().split('\t')
        bpos = int(bpos)
        freq = float(freq)
        if prot not in pwms:
            pwms[prot] = [[0.0]*4]
        elif bpos == len(pwms[prot]):
            pwms[prot].append([0.0]*4)
        pwms[prot][bpos][B2IND[base]] = freq
    fin.close()
    for p in pwms.keys():
        pwms[p] = np.array(pwms[p])
    return pwms

def makeLogos(pwms, logoDir):
    """
    Output function
    Place logos for every pwm and ouput to logoDir
    """
    keysToUse = sorted(pwms.keys())

    if not os.path.exists(logoDir):
        os.makedirs(logoDir)

    for k in keysToUse:
        pwm = pwms[k]
        logoLab = k
        makeNucMatFile('./tmp/','tmp',pwm)
        makeLogo('./tmp/tmp.txt',logoDir+k+'.pdf',
                 alpha = 'dna', colScheme = 'classic')
        os.system('rm %s' %'./tmp/tmp.txt')

def assignObsGrps(core, by = 'grpIDcore'):
    """
    Algorithm function.
    Assign observations to groups based on core similarity
    """
    grps = {}
    if by == 'none':
        # Each observation is assigned to its own group
        for k in core.keys():
            grps[k] = [k]
    elif by == 'grpIDcore':
        for k in core.keys():
            if not grps.has_key(core[k]):
                grps[core[k]] = [k]
            else:
                grps[core[k]].append(k)
    elif by == 'hd1':
        distMat, prots = getProtDistMat(core)

    return grps

def getFullX_fromFasta(fastaFile, aaPosList, edges):
    # Reads a fasta file of protein sequences and converts 
    # to the appropriate match state position representation,
    # and returns a corresponding X matrix

    # Convert to "core" contacting amino acids only
    corePos = [x - HMM_OFFSET for x in aaPosList]
    prot = readFromFasta(fastaFile)
    hmmerout, matchTab = './tmp/hasPWM.hmmer3.out.txt', './tmp/hasPWM.matchTab.txt'
    runhmmer3(hmmerout, HMM_FILE, fastaFile, HMM_NAME, getdescs(fastaFile), 
              hmmerDir = HMMER_HOME)
    core, full, trunc = \
        makeMatchStateTab(hmmerout, matchTab, prot, HMM_LEN, HMM_NAME, corePos = corePos)
    os.system('rm %s %s %s' %('./tmp/seq_fasta_hasPWM.fa', hmmerout, matchTab))
 
    # Remove examples for which a complete match couldn't be found
    keep = set([k for k in core.keys() if 'X' not in core[k] and '-' not in core[k]])
    subsetDict(core, keep)
    
    # Assign to distinct observation groups
    obsGrps = assignObsGrps(core, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots  

    nDoms = {}
    for p in uniqueProteins:
        nDoms[p] = len(core[p])/len(aaPosList)

    fullX, grpInd = formGLM_fullX(core, edges, uniqueProteins, obsGrps)
    return fullX, uniqueProteins, obsGrps, grpInd

def formGLM_fullX(core, edges, uniqueProteins, obsGrps, numAAs = 19, domainOrder = 1,
                  modelType = 'classifier'):
    """
    Function for forming a full X using all unique proteins including hold out proteins.
    :param core: a dictionary {"protein": "amino acid sequence"}
    :param edges: a dictionary {"dna base position": [amino acid positions]} defining the contact map
    :param uniqueProteins: an array of unique proteins
    :param domainOrder: set to -1 for arrayed domain proteins multidomain model is set up C-N w.r.t. binding site orientation, 1 otherwise 
    :return: a dictionary, each base position j corresponds to a matrix, 
    whose dimension is n*4 by E_j*numAAs, where n is the total number of domains across all proteins,
    and E_j the number of amino acids contact at the jth base position of the contact map (i.e., edges).
    """

    # Get the number of amino acid positions in the contact model
    allAApos = set()
    for j in range(len(edges)):
        for k in edges[j]:
            allAApos.add(k)
    coreLen = len(allAApos)

    X = {}
    grpInd = {}  # Maps obsGrp keys to (start, end) matrix row index pairs 
    for j in range(len(edges)):
        aa_pos = edges[j]
        E_j = len(aa_pos)
        X[j] = np.array([], dtype=np.int64).reshape(0, numAAs*E_j)
        rowNum = 0
        for g in obsGrps.keys():
            grpStart = rowNum
            for p in obsGrps[g]:
                nCores = len(core[p])/coreLen
                # Add X vectors for domains in reverse order for ZFs
                if domainOrder == -1:
                    cRange = range(nCores-1,-1,-1)
                else:
                    cRange = range(nCores)
                for c in cRange:
                    core_seq = core[p][coreLen*c:coreLen*(c+1)]
                    x = []
                    for i in aa_pos:
                        cur_x = [0] * numAAs
                        aa = A2IND[core_seq[i]]
                        if aa < numAAs:
                            cur_x[aa] = 1
                        x += cur_x
                    # For each protein, it has four lines of same X, 
                    # i.e, four samples for each protein in the logistic regression
                    if modelType == 'classifier':
                        x = np.tile(x, (4,1))
                        X[j] = np.concatenate((X[j], x), axis=0)
                        rowNum += 4             
            grpEnd = rowNum-1
            grpInd[g] = (grpStart, grpEnd)
    return X, grpInd

def formGLM_testX(X, startInd, endInd, modelType = 'classifier'):
    """
    Function for extracting X part from a full X just for hold out protein.
    :param X: a full X
    :param index: the index of the hold out domain(s) in the original X matrix
    :return: a dictionary. Each base position j corresponds to the X part for hold out protein.
    X[j] is a 4 by E_j*19 matrix, where E_j the number of amino acids contact to the jth base position.
    """

    if modelType == 'classifier':
        testX = {}
        for j in range(MWID):
            testX[j] = X[j][startInd:(endInd+1),]
    return testX  

def createGLMModel(trainX, trainY, trainW):
    """
    Function for performing eight weighted multinomial logistic regressions for eight base positions
    :param trainX: a dictionary of {base position: X matrix}
    :param trainY: a dictionary of {base position: y vector}
    :param trainW: a dictionary of {base position: weights, same length as y vector}
    :return:
    model: a dictionary of {base position: GLM model}
    """
    model = {}
    for j in range(MWID):
        #clf = LogisticRegression(fit_intercept=True, random_state=0, multi_class='multinomial', solver='newton-cg')
        clf = LogisticRegression(fit_intercept=True, random_state=0, 
                                 multi_class='multinomial', solver='newton-cg', 
                                 C = 1e9)
        model[j] = clf.fit(trainX[j], trainY[j], sample_weight=trainW[j])
    return model

def formGLM_Y(keysToUse, nDoms):
    """
    Function for forming Y for either train or test Y, each protein is [A,T,C,G], which is [0,1,2,3] in code
    :param keysToUse: an array of proteins used
    :return: an array of repeated [0,1,2,3]. And note that len(Y)/4 is the number of proteins used.
    """
    Y = {}
    for j in range(MWID):
        Y[j] = np.array([0,1,2,3] * np.sum([nDoms[k] for k in keysToUse]))
    return Y

def formGLM_trainW(pwms, uniqueProteins, nDoms, start, rev, modelType = 'classifier'):
    """
    Since we are using weighted multinomial logsitic regression, using only a single instance of each
    domain-DNA interface observation for each possible base outcome, this function creates weights for
    each row of the covariate matrix X based on given position weight matrix, starting position, 
    and orientation each PWM relative to the contact map.
    :param pwms: a dictionary of {"protein": position weight matrix}
    :param uniqueProteins: an array of unique proteins
    :param S: a dictionary {"protein": starting position}
    :param O: a dictionary {"protein": orientation \in {0,1}}
    :return: W: a dictionary. Each base position j corresponds to an array of weights,
    every four numbers represent one protein's weights.
    """

    # The PWM columns spanned by the (potentially multi-domain) interface under the current alignment
    weights = {}
    for i, protein in enumerate(uniqueProteins):
        if rev[protein] == 1:
            pwm = matrix_compl(pwms[protein])
        else:
            pwm = pwms[protein]
        weights[protein] = {}
        # Allows arrayed multi-domain proteins with overlaps
        for d in range(nDoms[protein]):
            for j in range(MWID):
                #print protein, start[protein], d, j, j+start[protein]+d*(MWID-RIGHT_OLAP), pwm[j+start[protein]+d*(MWID-RIGHT_OLAP)][:]
                if d == 0:
                    weights[protein][j] = pwm[j+start[protein]][:]
                else:
                    weights[protein][j] = \
                        np.concatenate((weights[protein][j],
                                        pwm[j+start[protein]+d*(MWID-RIGHT_OLAP)][:]), axis = None)

    W = {}
    if modelType == 'classifier':
        for j in range(MWID):
            W[j] = {}
            for i, protein in enumerate(uniqueProteins):
                W[j] = np.concatenate((W[j], weights[protein][j]), axis=None)
            W[j] = W[j][1:len(W[j])]
            # normalize weights for each W[j]
            W[j] = W[j] / sum(W[j])
            #print len(W[j])

    return W

def form_model(fullX, uniqueProteins, nDoms, pwms, start, rev):
    X = fullX
    Y = formGLM_Y(uniqueProteins, nDoms)
    W = formGLM_trainW(pwms, uniqueProteins, nDoms, start, rev)
    model = createGLMModel(X, Y, W)
    return model

def getPrecomputedInputs():
    ##############################################
    # Computes necessary data structures given PROT_SEQ_FILE, PWM_INPUT_TABLE,
    # HMM_FILE, HMM_LEN, HMM_NAME, HMM_OFFSET, and CONTACT_MAP.
    # See reference file formats in '../precomputedInputs/' (referenced at top of this script).
    
    # Get the PWM and protein info
    pwms = readPWMtab(PWM_INPUT_TABLE)
    prot = readFromFasta(PROT_SEQ_FILE)
    subsetDict(prot, set(pwms.keys()))
    
    # Get the contact map info
    edges_hmmPos = {}
    fin = open(CONTACT_MAP, 'r')
    fin.readline()
    aaPosList = set()
    for line in fin:
        bpos, aapos = tuple([int(x) for x in line.strip().split()])
        if bpos in edges_hmmPos:
            edges_hmmPos[bpos].append(aapos)
        else:
            edges_hmmPos[bpos] = [aapos]
        aaPosList.add(aapos)
    fin.close()
    aaPosList = sorted(list(aaPosList))
    edges = {}
    for bpos in edges_hmmPos:
        for aapos in edges_hmmPos[bpos]:
            if bpos in edges.keys():
                edges[bpos].append(aaPosList.index(aapos))
            else:
                edges[bpos] = [aaPosList.index(aapos)]

    # Convert to "core" contacting amino acids only
    corePos = [x - HMM_OFFSET for x in aaPosList]
    fasta = './tmp/seq_fasta_hasPWM.fa'
    writeToFasta(prot,'./tmp/seq_fasta_hasPWM.fa')
    hmmerout, matchTab = './tmp/hasPWM.hmmer3.out.txt', './tmp/hasPWM.matchTab.txt'
    runhmmer3(hmmerout, HMM_FILE, fasta, HMM_NAME, getdescs(fasta), 
              hmmerDir = HMMER_HOME)
    core, full, trunc = \
        makeMatchStateTab(hmmerout, matchTab, prot, HMM_LEN, HMM_NAME, corePos = corePos)

    # Removes test proteins/pwms from the dataset
    testProts = None
    if (EXCLUDE_TEST):
        fin = open(TEST_PROT_FILE,'r')
        testProts = [x.strip() for x in fin.readlines()]
        fin.close()
        for x in [pwms, core, full]:
            allProts = x.keys()
            for prot in testProts:
                try:
                    del x[prot]
                except KeyError:
                    next

    return pwms, core, full, edges, edges_hmmPos, aaPosList, testProts