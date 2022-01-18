# Finds an alignments of PWMs to that optimizes for dependency
# between amino acids occupying match state positions in an HMM
# and base positions in the set of PWMs.

from runhmmer import *
from matAlignLib import *
from pwm import makeNucMatFile, makeLogo
from copy import deepcopy
import numpy as np
import random
import multiprocessing
from getHomeoboxConstructs import getUniprobePWMs, getFlyFactorPWMs
from getHomeoboxConstructs import parseNoyes08Table, makeMatchStateTab
from getHomeoboxConstructs import readFromFasta, writeToFasta, subsetDict
import time
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

ORIENT = {'b08': {'Pitx1': 1}}   # Reference start position for one of the structures

OUTPUT_ALI = True  # Output properly aligned sequence logos for the best model
OUTPUT_MODEL = True  # Output the conditional and bg probs for best model
COMPARE_RUNS = True  # Output the model comparisons across multiple markov chains
MAKE_PREDS = True   # Make predicted PWMs for each protein based on best model
                     # and compare them to the input PWMs

# These are paramters for adjusting the model topology
# (i.e., the binarized contact matrix)
CORE_POS = 'useStructInfo'    # Uses the precomputed structural alignments 
OBS_GRPS = 'grpIDcore'        # Perform group updates based on common "core" AAs in proteins
APOS_CUT = 'cutAApos_1.0_0.05'  # Controls AA contact threshold for binarization
EDGE_CUT = 'edgeCut_1.0_0.05'   # Controls AA contact threshold for binarization
MAX_EDGES_PER_BASE = None       # None means to just ignore this parameter

TST_SET = 'b08'  # Just controls where code looks for the PCMs

MWID = 6  # Number of base positions in the PDSIM
N_CHAINS = 1 # Number of Markov chains if running in parallel
RAND_SEED = 382738372 #78374223 # Random seed for reproducibilty

# The set of "canonical" AA contacting positions, according to literature
#CANON9 = [2,3,5,6,47,50,51,54,55]

def getContactFractions(fname):
    """ Returns a nested dictionary containing the weighted
    fraction of observed structures for which each amino acid
    position was observed to contact at least one base.
    """
    fin = open(fname, 'r')
    fin.readline() # strip header
    wts = {'base': {}, 'backbone': {}}
    for line in fin:
        l = line.strip().split('\t')
        apos, cType, w = eval(l[0]), l[1], eval(l[2])
        wts[cType][apos] = w
    fin.close()

    return wts

def getEdgeWtMats(fname):
    """ Returns a matrix encoding the uniqueness weighted fraction
    of structures for which contact is made between each aa and base
    position (either backbone or base contact )
    """

    fin = open(fname, 'r')
    fin.readline() # strip header
    wts = {'base': {}, 'backbone': {}}
    allApos = set()
    allBpos = set()
    for line in fin:
        l = line.strip().split('\t')
        apos, bpos, w = tuple([eval(x) for x in l[:3]])
        cType = l[-1]
        if not (apos in wts[cType]):
            wts[cType][apos] = {}
        wts[cType][apos][bpos] = w
        allBpos.add(bpos)
        allApos.add(apos)
    fin.close()

    wtMats = {}
    for cType in ['base', 'backbone']:
        wtMats[cType] = \
            np.zeros((len(allApos),len(allBpos)), dtype = 'float')
        for i, apos in enumerate(sorted(allApos)):
            for j, bpos in enumerate(sorted(allBpos)):
                wtMats[cType][i,j] = wts[cType][apos][bpos]

    return wtMats, sorted(allApos), sorted(allBpos)

def getAAposByStructInfo(fname, wtCutBB, wtCutBase):
    """ Determine the amino acid positions to consider based on
    structural considerations.
    - fname specifies a text table of precomputed contact info
    - wtCutBB, and wtCutBase are the minimum weighed fraction of
      co-crystal structures in which an apos must contact DNA
      backbone or base, respectively
    """
    cWt = getContactFractions(fname)
    apos = list(set(cWt['backbone'].keys()) & set(cWt['base'].keys()))
    aaPosList = set([x for x in cWt['backbone'].keys() \
                     if cWt['backbone'][x] >= wtCutBB]) | \
                set([x for x in cWt['base'].keys() \
                     if cWt['base'][x] >= wtCutBase])
    return sorted(list(aaPosList))

def getEdgesByStructInfo(fname, aaPos, maxMwid, wtCutBB, wtCutBase,
                         maxEdgesPerBase = None, N51A_bpos = 2):
    """ Determine the bpos to apos dependency edges based on
    structural considerations.
    - fname specifies a text table of precomputed contact info
    - aaPos is a list (or set) of allowable apos
    - maxMwid is the maximum allowable motif width
    - wtCutBB, and wtCutBase are the minimum weighed fraction of
      co-crystal structures in which an apos x must contact
      bpos y on the backbone or base, respectively
    - N51A_bpos specifies the desired 0-indexed position for the
      adenine that interacts with N51 most strongly according
      (position 3 according to Christensen/Noyes, Cell 2008)
    - if maxEdgesPerBase is not None, it defines the maximum allowable 
      dependency edges per base positions, where the edges used are
      prioritized based on frequency the of base-contacting amino 
      acid positions (after removing aa-positions and edges based 
      on aaPos, wtCutBB, wtCutBase)
    """

    wtMats, rLabs, cLabs = getEdgeWtMats(fname)

    # Make N51 -> Ade contact be position 3
    mStart = np.argmax(wtMats['base'][rLabs.index(51),:]) - N51A_bpos
    edges_hmmPos = {}
    edgeWts = {}
    for bpos in range(maxMwid):
        edges_hmmPos[bpos] = []
        edgeWts[bpos] = []
    for apos in aaPos:
        for bpos in range(maxMwid):
            j = bpos+mStart
            bbWt = wtMats['backbone'][rLabs.index(apos),j]
            baseWt = wtMats['base'][rLabs.index(apos),j]
            if bbWt >= wtCutBB or baseWt >= wtCutBase:
                edges_hmmPos[bpos].append(apos)
                edgeWts[bpos].append(baseWt)

    # Restrict number of edges if desired
    if maxEdgesPerBase is not None:
        for bpos in edges_hmmPos.keys():
            keep = sorted(zip(edgeWts[bpos],edges_hmmPos[bpos]),
                          reverse = True)[:maxEdgesPerBase]
            edges_hmmPos[bpos] = sorted([x[1] for x in keep])

    # Ignores amino acids that do not contact any bases
    aaPosList = set()
    for k in edges_hmmPos.keys():
        aaPosList |= set(edges_hmmPos[k])
    aaPosList = sorted(list(aaPosList))


    # 0-index the edges with non-contacting amino acid positions removed
    edges = {}
    for bpos in range(maxMwid):
        edges[bpos] = []
        for x in edges_hmmPos[bpos]:
            edges[bpos].append(aaPosList.index(x))

    return edges, edges_hmmPos, aaPosList


def getHomeoboxData(dsets, maxMwid, aaPosList = [2,3,5,6,47,50,51,54,55],
                    aposCut = 'cutAApos_X', edgeCut = 'cutEdge_X',
                    maxEdgesPerBase = None, N51A_bpos = 2):
    """ Returns a set of amino acids occupying match states and
    mapped to the corresponding PWMs.

    aaPosList is a list of match states in the hmm to consider.
    if this value is set insead to 'useStructInfo', then the list
    is generated based on observed structural contacts
    """

    NOYES08_INFO = '../homeodomain/flyFactorSurvey/noyes_cell08/' + \
    '1-s2.0-S009286740800682X-mmc2_edited.csv'
    NOYES08_PWM_SRC = '../pwms/homeodomain/flyFactorSurvey_20180126.txt'
    BERGER08_FASTA = '../homeodomain/uniprobe_homeodomains_mouse.fa'
    BERGER08_PWM_SRC = '../pwms/homeodomain/uniprobe_20180126_Cell08/'
    HBOX_HMM = '../pfamHMMs/Homeobox.hmm'
    HBOX_HMM_LEN = 57
    HMM_NAME = 'Homeobox'
    HMM_OFFSET = 2

    # Optionally, get the core positions via structural analysis
    #print aposCut[-1]
    #print [eval(x) for x in aposCut.split('_')[1:]]
    if aaPosList == 'useStructInfo':
        assert aposCut[-1] != 'X'
        cutBB, cutBase = [eval(x) for x in aposCut.split('_')[1:]]
        wtFile = '../structuralAlignmentFiles/'+ \
            'Homeobox_weightedSum_distCut3.6_unionBases.txt'
        aaPosList = getAAposByStructInfo(wtFile, cutBB, cutBase)

    # Optionally, restrict edges based on structural info
    #print edgeCut[-1]
    #print [x for x in edgeCut.split('_')[1:]]
    if edgeCut[-1] != 'X':
        cutBB, cutBase = [eval(x) for x in edgeCut.split('_')[1:]]
        wtFile = '../structuralAlignmentFiles/' + \
            'Homeobox_weightedSum_distCut3.6.txt'
        edges, edges_hmmPos, aaPosList = \
            getEdgesByStructInfo(wtFile, aaPosList, maxMwid, cutBB,
                                 cutBase, maxEdgesPerBase = maxEdgesPerBase,
                                 N51A_bpos = N51A_bpos)
    else:
        # Includes all possible bpos-to-apos edges/dependencies
        edges = {}
        edges_hmmPos = {}
        for i in range(maxMwid):
            edges[i] = []
            edges_hmmPos[i] = []
            for j, apos in enumerate(aaPosList):
                edges[i].append(j)
                edges_hmmPos[i].append(apos)

    # Index relative to first match state of interest
    corePos = [x - HMM_OFFSET for x in aaPosList]

    seqs = {}   # The full length protein construct
    pwms = {}   # The full length pwm
    core = {}   # Concatenated "core" of extended match states for domain HMM
    full = {}   # Full extended match states for domain HMM
    trunc = {}  # Tells whether or not the HMM match was truncated
    for dset in dsets:

        # Get the sequences and PWMs for the particular datasets
        if dset == 'n08':
            seqs[dset] = parseNoyes08Table(NOYES08_INFO)
            pwms[dset] = getFlyFactorPWMs(NOYES08_PWM_SRC,
                                          subset = set(seqs[dset].keys()),
                                          whichPWM = 'Cell', countSmooth = 1)
            fstem = '/'.join(NOYES08_INFO.split('/')[:-1]) + \
                '/ffs_homeodomains_fly'
        elif dset == 'b08':
            seqs[dset] = readFromFasta(BERGER08_FASTA)
            pwms[dset] = getUniprobePWMs(BERGER08_PWM_SRC,
                                         subset = set(seqs[dset].keys()))
            print(pwms[dset])
            fstem = '/'.join(BERGER08_FASTA.split('/')[:-1]) + \
                '/uniprobe_homeodomains_mouse'


        # Subset seqeuences to those for which we have PWMs and
        # run hmmer on the set of sequences remaining
        fasta, hmmerout, matchTab = \
            fstem+'_hasPWM.fa', fstem+'_hasPWM.hmmer3.out.txt', \
            fstem+'_hasPWM.matchTab.txt'
        subsetDict(seqs[dset], set(pwms[dset].keys()))
        writeToFasta(seqs[dset], fasta)
        runhmmer3(hmmerout, HBOX_HMM, fasta, HMM_NAME, getdescs(fasta))
        core[dset], full[dset], trunc[dset] = \
            makeMatchStateTab(hmmerout, matchTab, seqs[dset],
                              HBOX_HMM_LEN, HMM_NAME, corePos = corePos)
        subsetDict(pwms[dset], set(core[dset].keys()))

    return seqs, pwms, core, full, trunc, aaPosList, edges, edges_hmmPos

def initStarts(pwms, mWid):
    """ Returns initial values for the starting positions and orientations
    of the pwms (random currently)
    """
    start, rev = {}, {}
    for k in pwms.keys():
        rev[k] = np.random.randint(low = 0, high = 2)
        start[k] = np.random.randint(low = 0,
                                     high = len(pwms[k])-mWid + 1)
    return start, rev

def getOrientedPWMs(pwms, rev, reorient = 1):
    """ Returns PWMs according to rev"""

    # Orient pwms according to rev
    pwms_o = {}
    for k in pwms.keys():
        if rev[k] == reorient:
            pwms_o[k] = matrix_compl(pwms[k])
        else:
            pwms_o[k] = deepcopy(pwms[k])
    return pwms_o

def getBackgroundAAprobs(seqs, smooth = 1, norm  = True,
                         logPr = False):
    """ Returns the per-position amino acid probabilities across
    the set of sequences using a smoothing count
    """
    p = {}
    for i in range(len(seqs[0])):
        p[i] = np.zeros(len(AMINO), dtype = 'float')
        for s in seqs:
            p[i][A2IND[s[i]]] += 1
        p[i] += smooth
        if norm:
            p[i] /= p[i].sum()
        if logPr:
            p[i] = np.log(p[i])
    return p

def getBackgroundBaseProbs(pwms, starts, mWid, smooth = 0.01,
                           norm = True, logPr = False):
    """ Returns the per-position base probabilities across
    the set of pwms
    """
    p = {}
    for i in range(mWid):
        p[i] = np.zeros(len(BASE), dtype = 'float')
        for j, x in enumerate(pwms):
            p[i] += x[starts[j]+i,:]
        p[i] += smooth
        if norm:
            p[i] /= p[i].sum()
        if logPr:
            p[i] = np.log(p[i])
    return p

def getCondModel(pwms, seqs, start, edges, keysToUse = None,
                  smooth = 0.01, logProbs = False):
    """ Returns a smoothed estimate of the cond P(a_i|b_j) model,
    given the start positions and pwm orientations.

    - pwms, seqs, and start are dicts indexed by interface IDs
    - edges maps base positions to dependent aa positions
    - if logProbs then return conditional probs as logProbs
    """

    # For for complete bipartite edge set
    """
    edges = {} # bpos -> abpos
    for i in range(mWid):
        edges[i] = []
        for j in range(nApos):
            edges[i].append(j)
    """

    # Compute per-position background probs
    aaProbs = getBackgroundAAprobs([seqs[k] for k in keysToUse],
                                   logPr = False)

    nApos = len(seqs[seqs.keys()[0]])
    mWid = len(edges.keys())

    # Maps each amino-base edge to a 4X20 conditional
    # prob distributions based on pwms
    p = {}
    p_ba = {}  # Estimate for P(b_j|a_i)
    for j in range(mWid):
        for i in edges[j]:
            p[i,j] = np.zeros((len(AMINO), len(BASE)),
                              dtype = 'float')
            p_ba[i,j] = np.zeros((len(AMINO), len(BASE)),
                                 dtype = 'float')

            # Compute P(b_j|a_i) as smoothed average across PWMs
            for k in keysToUse:
                aa = seqs[k][i]
                p_ba[i,j][A2IND[aa],:] += pwms[k][start[k]+j,:]
            for aa in A2IND.keys():
                p_ba[i,j][A2IND[aa],:] += smooth
                p_ba[i,j][A2IND[aa],:] /= p_ba[i,j][A2IND[aa],:].sum()

            # Use Bayes' rule to convert to P(a_i|b_j)
            for aa in A2IND.keys():
                p[i,j][A2IND[aa],:] = \
                    p_ba[i,j][A2IND[aa],:]*aaProbs[i][A2IND[aa]]
            for base in B2IND.keys():
                p[i,j][:,B2IND[base]] /= p[i,j][:,B2IND[base]].sum()

            # Switch to log-probs if desired
            if logProbs:
                p[i,j] = np.log(p[i,j])

    # Compute a marginal P(b_j) for each j
    bProbs = {}
    for j in range(mWid):
        bProbs[j] = np.zeros(len(BASE), dtype = 'float')
        for i in edges[j]:
            tmp = np.zeros(len(BASE), dtype = 'float')
            for aa in A2IND.keys():
                tmp += p_ba[i,j][A2IND[aa],:]*aaProbs[i][A2IND[aa]]
            bProbs[j] += tmp/tmp.sum()
        bProbs[j] /= bProbs[j].sum()

    # Switch to log-probs if desired
    if logProbs:
        for k in aaProbs.keys():
            aaProbs[k] = np.log(aaProbs[k])
        for k in bProbs.keys():
            bProbs[k] = np.log(bProbs[k])

    return p, aaProbs, bProbs

def getLogLikelihood(pwm, aSeq, cond, marg, s, edges):
    """ Returns a value proportional (in the limit) to the conditional
    log-likelihood of the mWid pwm and starting at position s, given
    the amino acid sequence aSeq, current estimated conds P(a_i|b_j),
    and current estimated marg P(b_j).
    The likelihood is computed under the assumptions that:
    1. - all base positions are independent of one another, and
    2. - all amino acid positions are conditionally independent given bases.

    I.e, the complete conditional likelihood is the product of likelihoods
    of Naive Bayes' models for P(b_j|a_i) across all i.
    """

    mWid = len(edges.keys())
    sampSz = 100
    loglik = 0.0
    for j in range(mWid):
        for b in B2IND.keys():
            tmp = marg[j][B2IND[b]]
            for i in edges[j]:
                tmp += cond[i,j][A2IND[aSeq[i]],B2IND[b]]
            tmp *= int(round(sampSz*pwm[j+s][B2IND[b]]))
            loglik += tmp
    return loglik

def getLLsum(pwms, aSeqs, cond, bProbs, start, edges, keysToUse = None):
    """
    Compute the summed log likelihoods across all the data
    """

    if keysToUse is None:
        keysToUse = pwms.keys()
    ll = 0.0
    for k in keysToUse:
        ll += getLogLikelihood(pwms[k], aSeqs[k], cond, bProbs,
                               start[k], edges)
    return ll

def sampleStartPos(pwm, aSeq, cond, marg, edges):
    """ Returns a sampled new start position and orientation for pwm
    given aSeq, the cond P(a_i|b_j), and the marg P(b_j).
    """

    # Construct the sampling distirbution
    seqWid = len(pwm)
    mWid = len(edges.keys())
    lls = []
    for pwm_o in [pwm, matrix_compl(pwm)]:
        for s in range(seqWid-mWid):
            ll = getLogLikelihood(pwm_o, aSeq, cond, marg, s, edges)
            lls.append(ll)

    #print lls
    sdistr = np.exp(lls - max(lls)) # Standardize to most likely start site
    sdistr = sdistr/sdistr.sum()

    # Sample and return the new start position/orientation
    cdistr = np.cumsum(sdistr)
    rand = np.random.random_sample()
    i = 0
    while i < len(cdistr):
        if rand <= cdistr[i]:
            break
        i += 1
    s = i

    rev = 0
    if s >= len(lls)/2:
        rev = 1
    s = s%(len(lls)/2)

    return s, rev, lls

def reverseBaseOrient(pwms, mWid = None, start = None, rev = None):
    """ Reverses the orientation of all pwms and optionally also
    corrects starting positions and reverse flags if a motif width
    is provided.
    """
    pwms_o = {k: matrix_compl(pwms[k]) for k in pwms.keys()}
    if mWid is not None:
        opp = {0:1,1:0}
        start_o = {k: len(pwms[k])-start[k]-mWid for k in start.keys()}
        rev_o = {k: opp[rev[k]] for k in rev.keys()}
        return pwms_o, start_o, rev_o
    else:
        return pwms_o

def gibbsSample(pwms, aSeqs, edges, obsGrps, verbose = False,
                maxIters = 25,randSeed = None,orientKey = None,orient = 1):
    """ Runs a Gibbs sampler until convergence and return the
    the resulting model info as a tuple
    """

    # Set the random seed again in case this is called by a parent process
    if randSeed is not None:
        np.random.seed(randSeed)
    else:
        np.random.seed()

    mWid = len(edges.keys())
    # Init the latent variables randomly
    start, rev = initStarts(pwms, mWid)

    # Alternate between model contruction and latent variable updates
    # until the latent variables cease to change
    converged = False
    nIters = 0
    if verbose:
        print("nIters: %d" %nIters)
        print("\tstarts:", start.values()[:20])
        print("\trevers:", rev.values()[:20])


    dir = "debug_res_naiveBayes"
    if not os.path.exists(dir):
        os.makedirs(dir)


    while not converged and nIters < maxIters:

        valuesChanged = 0
        # A single iteration, allowing group holdouts
        for grpKey in np.random.permutation(obsGrps.keys()):
            grp = obsGrps[grpKey]
            keysToUse = set(pwms.keys())-set(grp)

            # Update the conditional probability (A_i|B_j) tables,
            # excluding the held out group of observations
            cond, aaProbs, bProbs = \
                getCondModel(getOrientedPWMs(pwms, rev),aSeqs,start,edges,
                             keysToUse=keysToUse,logProbs=True)

            # Sample a new starting position for each held out interface
            for holdout in np.random.permutation(grp):
                s, r, lls = sampleStartPos(pwms[holdout],aSeqs[holdout],cond,
                                      bProbs,edges)
                if s != start[holdout] or r != rev[holdout]:
                    valuesChanged += 1
                start[holdout] = s
                rev[holdout] = r

        nIters += 1

        if verbose:
            print("nIters: %d" %nIters)
            print("\tstarts:", start.values()[:20])
            print("\trevers:", rev.values()[:20])
            print("\tvaluesChanged: ", valuesChanged)
        # Check for convergence
        if not valuesChanged:
            converged = True

    # Compute the logLikelihood of all the data under all the
    # (converged) hidden params
    cond, aaProbs, bProbs = \
        getCondModel(getOrientedPWMs(pwms, rev), aSeqs, start, edges,
                     keysToUse = pwms.keys(),logProbs = True)
    ll = getLLsum(getOrientedPWMs(pwms, rev), aSeqs, cond, bProbs,
                  start, edges)

    # Reverse all orientations if the key pwm is reversed
    pwms_o = getOrientedPWMs(pwms, rev)
    start_o, rev_o = deepcopy(start), deepcopy(rev)
    reorient = False
    if rev[orientKey] != orient:
        pwms_o, start_o, rev_o = \
            reverseBaseOrient(pwms_o, mWid = mWid, start = start,rev = rev)
        reorient = True

    # Recompute probs without logging to return
    cond, aaProbs, bProbs = \
        getCondModel(pwms_o, aSeqs, start_o, edges, keysToUse = pwms_o.keys())

    return {'aaProbs': aaProbs, 'bProbs': bProbs, 'cond': cond,
        'start': start, 'rev': rev, 'll': ll, 'nIter': nIters,
        'seed': randSeed, 'reorient': reorient}

def runGibbs(pwms, core, edges, obsGrps, verbose = False,
             kSamps = 25, orientKey = None, orient = None):
    """ Runs the gibbsSampler routine K times using parallel executions
    """

    ncpus = multiprocessing.cpu_count()-1
    maxiter = 25
    p = multiprocessing.Pool(ncpus)
    procs = []
    for k in range(kSamps):
        args = (pwms, core, edges, obsGrps, verbose,
                maxiter,np.random.randint(0,1e9),
                orientKey, orient,)
        procs.append(p.apply_async(gibbsSample, args=args))
    res = [x.get() for x in procs]
    return res

#### Output/evaluation functions ###

def getAlignedPWMs(pwms, aSeqs, start, mWid, flipAli = False):
    """ Returns a new set of PWMs, truncated on each side to the
    aligned region
    """

    npwms = {}
    for p in pwms.keys():
        pwm = pwms[p]
        npwm = np.zeros((mWid,4), dtype = 'float')
        s = start[p]
        for i in range(mWid):
            if i+s > len(pwm) - 1:
                npwm[i,:] = 4*[0.25]
            else:
                npwm[i,:] = pwm[i+s,:]
        if flipAli:
            npwm = matrix_compl(npwm)
        npwms[p] = npwm
    return npwms

def makeAllLogos(pwms, aSeqs, logoDir, keysToUse = None):
    """ Place logos for every pwm and ouput to logoDir
    """

    if keysToUse == None:
        keysToUse = sorted(pwms.keys())

    for k in keysToUse:
        pwm, core = pwms[k], aSeqs[k]
        logoLab = '_'.join([k, core])
        makeNucMatFile('./tmp/','tmp',pwm)
        makeLogo('./tmp/tmp.txt',logoDir+logoLab+'.pdf',
                 alpha = 'dna', colScheme = 'classic')
        os.system('rm %s' %'./tmp/tmp.txt')

def outputRunSummary(fname, start, rev, nIter, ll, pwmWid,
                     mWid, seed, reorient, core):
    """ Output information to compare the K runs of the Gibbs sampler
    """

    fout = open(fname, 'w')
    fout.write('nRun\tprot\tstart\trev\tnIter\tloglik\tpwmWid' + \
               '\tmWid\trSeed\treorient\tcore\n')
    for i in range(len(start)):
        for k in sorted(start[i].keys()):
            fout.write('%d\t%s\t%d\t%d\t%d\t%e\t%d\t%d\t%d\t%d\t%s\n' \
                       %(i,k,start[i][k],rev[i][k],nIter[i],
                         ll[i], pwmWid[k], mWid, seed[i],
                         int(reorient[i]), core[k]))
    fout.close()

def makeCondProbsTab(cond, fname):
    """ Create a table for summarizing the per-position-pair
    conditional probabilities P(a,b)_(i,j)
    """

    fout = open(fname, 'w')
    fout.write('apos\tbpos\taa\tbase\tprob\n')
    for k in sorted(cond.keys()):
        apos, bpos = k[0], k[1]
        x = cond[k]
        for aa in AMINO:
            for b in BASE:
                fout.write('%s\t%s\t%s\t%s\t%e\n' \
                           %(apos,bpos,aa,b,x[A2IND[aa],B2IND[b]]))
    fout.close()

def makeBackgroudAAprobTab(p, fname):
    fout = open(fname, 'w')
    fout.write('apos\taa\tprob\n')
    for i in sorted(p.keys()):
        for j in sorted(IND2A.keys()):
            aa = IND2A[j]
            prob = p[i][j]
            fout.write('%d\t%s\t%e\n' %(i, aa, prob))
    fout.close()

def makeBackgroudBaseProbTab(p, fname):
    fout = open(fname, 'w')
    fout.write('bpos\tbase\tprob\n')
    for i in range(len(p)):
        for j in sorted(IND2B.keys()):
            base = IND2B[j]
            fout.write('%d\t%s\t%e\n' %(i, base, p[i][j]))
    fout.close()

def makePWMtab(pwms, fname):
    fout = open(fname, 'w')
    fout.write('prot\tbpos\tbase\tprob\n')
    for k in pwms.keys():
        m = pwms[k]
        for i in range(len(m)):
            for j in range(len(m[i])):
                fout.write('%s\t%d\t%s\t%e\n' %(k,i,IND2B[j],m[i,j]))
    fout.close()

def getProtDistMat(core):
    """ Create a matrix of hamming distances between protein core seqs
    """

    prots = sorted(core.keys())
    protDist = np.zeros((len(prots),len(prots)), dtype = 'float')
    for i, k1 in enumerate(prots):
        c1 = core[k1]
        for j, k2 in enumerate(prots):
            c2 = core[k2]
            hdist = 0
            for xpos in range(len(c1)):
                if c1[xpos] != c2[xpos]:
                    hdist += 1
            protDist[i,j] = hdist
    return protDist, prots

def assignObsGrps(core, by = 'grpIDcore'):
    """ Assign observations to groups based on core similarity
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

def getClosestNoyesProt(c1, ref):
    """ Returns a dictionary mapping labels from dict c1 to
    a (label, distance) pair corresponding to an entry from
    the ref dict that is closest in terms of hamming distance
    """

    closest = {}
    for k in c1.keys():
        minDist = 1e9
        minKey = ''
        x = c1[k]
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



def main():
    # Reproducible randomness
    mWid = MWID
    np.random.seed(RAND_SEED)

    # Get the data
    if CORE_POS == 'canon9':
        aaPosList = CANON9
    else:
        aaPosList = CORE_POS

    # This sets reads in the PWMs and sets up the contact structure 
    # between DNA and protein (i.e., the edge set in the model).    
    seqs, pwms, core, full, trunc, aaPosList, edges, edges_hmmPos = \
        getHomeoboxData(['n08', 'b08'], MWID, aaPosList = aaPosList,
                        aposCut = APOS_CUT, edgeCut = EDGE_CUT,
                        maxEdgesPerBase = MAX_EDGES_PER_BASE,
                        N51A_bpos = (MWID-6)/2+2)

    ### NOTE!! - the code currently assumes all available pwms
    ### are at least as wide as mWid ... this needs to be fixed eventually

    # Test using just one PWM dataset at a time for now, eventually may combine dsets
    tstSet = TST_SET

    # Create the output directories
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
    mainOutDir = '../homeodomain/gibbsAlign_output%s' %(suffix)
    if not os.path.exists(mainOutDir):
            os.makedirs(mainOutDir)

    # Subset to the test group and get orientation of key PFM
    orientKey = ORIENT[tstSet].keys()[0]
    orient = ORIENT[tstSet][orientKey]
    seqs, pwms, core, full, trunc = seqs[tstSet], pwms[tstSet], \
        core[tstSet], full[tstSet], trunc[tstSet]
    keep = set([k for k in core.keys() if 'X' not in core[k]])
    [subsetDict(x, keep) for x in [seqs, pwms, core, full, trunc]]

    # Assign cores to similarity groups
    obsGrps = assignObsGrps(core, by = OBS_GRPS)
    with open(mainOutDir+'obsGrpTab.txt', 'w') as fout:
        for k in sorted(obsGrps.keys()):
            for x in obsGrps[k]:
                fout.write('%s\t%s\n' %(k,x))

    print("Output written to: %s" %mainOutDir)
    # Write the PWMs used and other basic info to files
    makePWMtab(pwms, mainOutDir+'pwmTab.txt')
    with open(mainOutDir+'coreSeqTab.txt', 'w') as fout:
        fout.write('\n'.join(['%s\t%s' %(k,core[k]) \
                             for k in sorted(core.keys())]))
    with open(mainOutDir+'aaPosList.txt', 'w') as fout:
        fout.write('\n'.join(['%d' %(k) for k in aaPosList]))
    with open(mainOutDir+'edgeList.txt', 'w') as fout:
        fout.write('\n'.join(['%d\t%s' %(k, str(edges_hmmPos[k])) \
                             for k in sorted(edges_hmmPos.keys())]))

    # How many times is each amino acids observed in each position?
    aaCounts = getBackgroundAAprobs(core.values(), smooth = 0, norm = False)
    #print aaCounts
    makeBackgroudAAprobTab(aaCounts, mainOutDir+'aaCounts.txt')

    #print edges
    #print obsGrps
    #print orientKey, orient
    #print core

    # For debugging/testing purposes use this - runGibbs calls gibbsSample 
    # function multiple times in parallel with the different burn in criteria
    #gibbsSample(pwms,core,edges,obsGrps, verbose = False,
    #            maxIters = 1,randSeed = np.random.randint(0,1e9),
    #            orientKey = orientKey, orient = orient)


    # Code used for the analysis.  I.e., runs multiple markov chains in parallel
    # with different random starting posiitons and outputs information for the 
    # one that obtained highest model likelihood.

    startTime = time.time()
    print("Running %d markov chains ..." %N_CHAINS)
    # Run the Gibbs sampler K times and compare the outputs
    res = runGibbs(pwms,core,edges,obsGrps, verbose = False,
                   kSamps = N_CHAINS, orientKey = orientKey, orient = orient)
    print("Ran in %.2f seconds" %(time.time()-startTime))
    aaProbs = [x['aaProbs'] for x in res]
    bProbs = [x['bProbs'] for x in res]
    cond = [x['cond'] for x in res]
    start = [x['start'] for x in res]
    rev = [x['rev'] for x in res]
    ll = [x['ll'] for x in res]
    nIter = [x['nIter'] for x in res]
    seed = [x['seed'] for x in res]
    reorient = [x['reorient'] for x in res]
    opt = np.argmax(ll)
    print("Max iter:", opt, max(ll))

    with open('nbres.pickle', 'wb') as f:
        pickle.dump(res, f)

    # Output aligned logos based on the optimal PWM model
    # for the top 2 most highest likelihood models

    # Output the final per-position(-pair) cond and marginal probability models
    if OUTPUT_MODEL:
        if not os.path.exists(mainOutDir):
            os.makedirs(mainOutDir)
        makeCondProbsTab(cond[opt], mainOutDir+'condProbs.txt')
        makeBackgroudBaseProbTab(bProbs[opt], mainOutDir+'baseProbs.txt')
        makeBackgroudAAprobTab(aaProbs[opt], mainOutDir+'aaProbs.txt')

    if COMPARE_RUNS:
        # Output the final conditional per-position conditional probability models
        cmpDir = mainOutDir+'compareRuns/'
        if not os.path.exists(cmpDir):
            os.makedirs(cmpDir)
        pwmWid = {k: len(pwms[k]) for k in pwms.keys()}
        outputRunSummary(cmpDir+'cmpRuns.txt', start, rev, nIter, ll,
                         pwmWid, mWid, seed, reorient, core)
        for i in range(len(res)):
            makeCondProbsTab(cond[i],cmpDir+'condProbs_%d.txt'%i)
            makeBackgroudBaseProbTab(bProbs[i], cmpDir+'baseProbs_%d.txt'%i)
            makeBackgroudAAprobTab(aaProbs[i], cmpDir+'aaProbs_%d.txt'%i)

    if OUTPUT_ALI:
        print("Creating aligned logos ...")
        for i, o in enumerate([opt]):
            logoDir = mainOutDir + 'logos_%d/' %(i+1)
            if not os.path.exists(logoDir):
                os.makedirs(logoDir)
            flipAli = False
            if reorient[o]:
                flipAli = True
            aliPWMs = getAlignedPWMs(getOrientedPWMs(pwms, rev[o]),
                                     core, start[o], mWid, flipAli = flipAli)
            print("logoDir", logoDir)
            makeAllLogos(aliPWMs, core, logoDir)

    if MAKE_PREDS:
        os.system('/System/Library/Frameworks/Python.framework/Versions/2.7/bin/python2.7  naiveBayes_predict.py %s' \
                  %(' '.join([TST_SET,CORE_POS,OBS_GRPS,
                             APOS_CUT, EDGE_CUT, str(MAX_EDGES_PER_BASE),
                             str(MWID),str(N_CHAINS),str(RAND_SEED)])))

if __name__ == '__main__':
    main()
