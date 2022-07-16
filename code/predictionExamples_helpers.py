# Helper functions and global variables for predictionExample.py

from sklearn.linear_model import LogisticRegression
from gibbsAlign_GLM import makePWMtab
from getHomeoboxConstructs import readFromFasta, writeToFasta, subsetDict, makeMatchStateTab
from pwm import makeNucMatFile, makeLogo, rescalePWM
from matAlignLib import matrix_compl
from runhmmer import runhmmer3, getdescs
import numpy as np
import pickle, os

DOMAIN_TYPE = 'homeodomain'  # Set to either 'zf-C2H2' or 'homeodomain'
OBS_GRPS = 'grpIDcore'

if DOMAIN_TYPE == 'zf-C2H2':
    PROT_SEQ_FILE = '../precomputedInputs/zf-C2H2/prot_seq_fewZFs_hmmerOut_clusteredOnly_removeTooShort.txt'  # Input protein domain file subsetted to relvant amino acid contacting positions
    PROT_SEQ_FILE_FFS = '../flyFactorSurvey/enuameh/enuameh_perFinger_processedProtInfo.txt'
    PWM_INPUT_TABLE = '../precomputedInputs/zf-C2H2/pwmTab_fewZFs_clusteredOnly_removeTooShort.txt'   # A table of PWMs corresponding to prots in PROT_SEQ_FILE
    PWM_INPUT_FILE_FFS = '../flyFactorSurvey/enuameh/flyfactor_dataset_A.txt'
    CONTACT_MAP = '../precomputedInputs/zf-C2H2/zf-C2H2_contactMap.txt'  # A contact map for the domain family
    SEED_FILE = '../flyFactorSurvey/enuameh/enuameh_startPosInfo.txt'    # Initial seeds based on Enuameh et al. 2013
    MWID = 4
    RIGHT_OLAP = 1     # Number of 3' bases in contact map overlapping with previous domain instance (if multi-domain) - 1 for zf-C2H2
    DOMAIN_ORDER = -1  # Since zf-C2H2s bind 3' to 5', consider them in reverse order
elif DOMAIN_TYPE == 'homeodomain':
    PROT_SEQ_FILE = '../precomputedInputs/proteins_homeodomains_hasPWM.fa'  # Input protein sequence fasta file
    PWM_INPUT_TABLE = '../precomputedInputs/pwmTab_homeodomains_all.txt'   # A table of PWMs corresponding to prots in PROT_SEQ_FILE
    CONTACT_MAP = '../precomputedInputs/homeodomain_contactMap.txt'  # A contact map for the domain family
    MWID = 6
    HMM_FILE = '../pfamHMMs/Homeobox.hmm'    # Location of hmm file
    HMM_LEN = 57                             # Number of match states in HMM file
    HMM_NAME = 'Homeobox'                    # A name for the HMM
    HMM_OFFSET = 2                           # Used to offset HMM states to canonical numbering scheme for Homoedomains
    HMMER_HOME = None
    EXCLUDE_TEST = False
    RIGHT_OLAP = 0     # No domain base overlap for single-domain proteins
    DOMAIN_ORDER = 1   # Standard domain ordering since only 1 domain per protein

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

def getFullX_fromFasta_HMMer3(fastaFile, aaPosList, edges):
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

def getFullX_fromFile_zfC2H2(infile, aaPosList, edges):
    # Read input zfC2H2s from special format file

    fin = open(infile, 'r')
    core = {}  # TO BE COMPLETED
    for line in fin:
        l = line.strip().split()
        prot, cs = l[0], l[1]
        if core.has_key(prot):
            core[prot] += cs
        else:
            core[prot] = cs

    # Assign to distinct observation groups
    obsGrps = assignObsGrps(core, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots  

    nDoms = {}
    for p in uniqueProteins:
        nDoms[p] = len(core[p])/len(aaPosList)
        #print p, core[p], nDoms[p]
    fullX, grpInd = formGLM_fullX(core, edges, uniqueProteins, 
                                  obsGrps, domainOrder = DOMAIN_ORDER)
    
    return fullX, uniqueProteins, obsGrps, grpInd, nDoms

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

def getPrecomputedInputs_zfC2H2(rescalePWMs = False, ffsOnly = False, includeB1H = False):
    # Get the necesarry precomputed information for the C2H2-ZF inputs
    # Assumes an input table mapping protein IDs to ordered arrays of domains
    # subsetted to the appropriate base-contacting positions according to HMMER v2.3.2.
    # ffsOnly using only the flyfactorsurvey data for debugging purposes

    def getPWM(fpath, tfs, motifs, verbose = False, ID_field = "TF Name"):
        # Extracts motifs from the CIS-BP PWM.txt file subset to 
        # only those whose TF_Name is in tfs (set) and whose 
        # Motif_ID is in motifs (set)
        fin = open(fpath)
        pwms = {}
        line = fin.readline()

        while line != "":
            lineArr = line.split("\t")
            if verbose:
                print lineArr
            if lineArr[0] == ID_field:
                tf = lineArr[1].rstrip()
            if lineArr[0] == "Pos":
                pwm = []
            if lineArr[0] == "Motif":
                motif = lineArr[1].rstrip()
            if len(lineArr) == 5 and lineArr[0] != "Pos":
                lineArr = np.array(lineArr)
                onevec_list = lineArr.astype(np.float)[1:5].tolist()
                pwm.append(onevec_list)
            if lineArr[0] == '\n' and tf in tfs and motif in motifs:
                pwms[tf] = np.array(pwm)
                line = fin.readline()
            line = fin.readline()
        return pwms

    def getProteinInfo_zfC2H2_FFS(infile):
        core = {}
        fin = open(infile, 'r')
        fin.readline()
        for line in fin:
            l = line.strip().split()
            pname, coreSeq = l[0], l[3]
            if core.has_key(pname):
                core[pname] += coreSeq
            else:
                core[pname] = coreSeq
        fin.close()    
        return core

    def getFlyFactorPWMs_zfC2H2(ffsPWMfile, prots = set(), smooth = 0):
        # Read the fly factor survey pwms into numpy arrays

        fin = open(ffsPWMfile, 'r')
        line = fin.readline()
        pwms = {}
        while line != '':
            if line[0] == '>':
                l = line.strip().split('\t')
                motif = l[0][1:]
                line = fin.readline()
                if motif.split('_')[0] in prots:
                    prot = motif.split('_')[0]
                else:
                    while line != '' and line[0] != '>':
                        line = fin.readline()
                    continue
                if prot+'_SOLEXA' in motif:
                    pwm = []
                    while line != '' and line[0] != '>':
                        pwm.append([float(x)+smooth for x in line.strip().split('\t')])
                        line = fin.readline()
                    pwms[prot] = np.array(pwm)
                else:
                    while line != '' and line[0] != '>':
                        line = fin.readline()
                    continue
        fin.close()

        # Normalize the motifs
        for p in pwms.keys():
            for i in range(len(pwms[p])):
                pwms[p][i] = pwms[p][i]/pwms[p][i].sum()    
        return pwms

    def readSeedAlignment(infile, include = set()):
        # Reads table of seed start/orientations (as created by getFixedStarts_fromStructures())
        fixedStarts = {}
        fin = open(infile)
        fin.readline()
        for line in fin:
            l = line.strip().split()
            if l[0] in include:
                fixedStarts[l[0]] = {'start': int(l[1]), 'rev': int(l[2])}
        fin.close()
        return fixedStarts

    # Get the PWM info
    # Get the motif IDs of interest
    fin = open('../cis_bp/C2H2-ZF/motifTable_mostRecent_fewZFs_clusteredOnly_removeTooShort.txt', 'r')
    pwms = {}
    if not ffsOnly:
        fin.readline()
        motifMap = {}
        motifs = set()
        for line in fin:
            l = line.strip().split('\t')
            prot, motif = l[0],l[3]
            motifMap[prot] = motif
            motifs.add(motif)
        fin.close()
        #print motifMap
        pwms = getPWM('../cis_bp/C2H2-ZF/PWM.txt', set(motifMap.keys()), motifs, ID_field = 'TF')
    
    # Read in the fly-factor survey info
    core_ffs = getProteinInfo_zfC2H2_FFS(PROT_SEQ_FILE_FFS)
    pwms_ffs = getFlyFactorPWMs_zfC2H2(PWM_INPUT_FILE_FFS, prots=set(core_ffs.keys()), smooth = 1)
    #print(len(core_ffs), len(pwms_ffs))
    subsetDict(core_ffs, set(pwms_ffs.keys()))
    subsetDict(pwms_ffs, set(pwms_ffs.keys()))
    #print(len(core_ffs), len(pwms_ffs))

    # Get protein info
    core = {}
    if not ffsOnly:
        fin = open(PROT_SEQ_FILE, 'r')
        fin.readline()
        for line in fin:
            l = line.strip().split()
            pname, coreSeq = l[0], l[7]
            if core.has_key(pname):
                core[pname] += coreSeq
            else:
                core[pname] = coreSeq
        fin.close()    
        subsetDict(core, set(pwms.keys()))

    # Combine to a single set of prots/PWMs
    for k in core_ffs.keys():
        core[k] = core_ffs[k]
        pwms[k] = pwms_ffs[k]

    # Read in the fly-factor survey info
    if includeB1H:
        core_b1h, pwms_b1h = getB1Hmotifs('../cis_bp/C2H2-ZF/B1H.motifs.long.pfm.scaled.noReps.txt')
        #print core_b1h
        # Combine to a single set of prots/PWMs
        for k in core_b1h.keys():
            core[k] = core_b1h[k]
            pwms[k] = pwms_b1h[k]
   
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

    # Rescale the PWMs?
    if rescalePWMs:
        for k in pwms.keys():
            pwms[k] = rescalePWM(pwms[k], maxBaseSelect = 50)

    return pwms, core, edges, edges_hmmPos, aaPosList

def predictSpecificity_array_ZF(fullX, model, startInd, arrayLen, wtB1 = 0.5,
                                rescaleIC = True):
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
    #print pwms
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

    if rescaleIC:
        return rescalePWM(np.array(pwm), maxBaseSelect = 50)
    else:
        return np.array(pwm)