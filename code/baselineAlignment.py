# Parse the output of STAMP to visualize aligned 6bp logos for those
# that match with exact core sequence to a protein from Noyes 08 with 
# known alignment.

import numpy as np
import os, re
from gibbsAlign_GLM import getHomeoboxData, makeAllLogos
from gibbsAlign_naiveBayes import getAlignedPWMs, getOrientedPWMs, reverseBaseOrient
from getHomeoboxConstructs import subsetDict
from GLM_predict import getNoyesSpecGrps, getClosestNoyesProt, bestAliScoreWithNoyes, NOYES_SPEC_GRPS

CORE_POS = 'useStructInfo'    # Uses the precomputed structural alignments
OBS_GRPS = 'grpIDcore'        # Perform group updates based on common "core" AAs in proteins
APOS_CUT = 'cutAApos_1.0_0.05'  # Controls AA contact threshold for binarization
EDGE_CUT = 'edgeCut_1.0_0.05'   # Controls AA contact threshold for binarization
MAX_EDGES_PER_BASE = None       # None means to just ignore this parameter
MWID = 6

STAMP_ALI_FILE = '../STAMP/Stamp_Results_Homeodomain_10-21-20_alignment.txt'
ORIG_CONSENSUS_FILE = '../STAMP/stampConsensusCodeLabels.txt'
COMP_BASE = {'A':'T','C':'G','G':'C','T':'A','M':'K','K':'M','R':'Y','W':'W','S':'S','Y':'R','-':'-','N':'N'}

def makeRevCompl(s):
    # Makes the reverse complement of an alignment string 
    # or a consensus sequence

    return ''.join([COMP_BASE[x] for x in s[::-1]])

def readSTAMPaliFile(stampAlignmentFile, rc_flag = False):
    # Returns a dictionary mapping each protein to its alignment
    # string from the STAMP alignemnt

    fin = open(stampAlignmentFile, 'r')
    aliString = {}
    #consensusStr = {}
    for line in fin:
        l = line.strip().split()
        prot, ali = l[0][:-1], l[1]
        if rc_flag:
            ali = makeRevCompl(ali)
        aliString[prot] = ali
        #consensusStr[prot] = ali.lstrip('-').strip('-')
    fin.close()
    return aliString

def readForwardConsensus(origConsensusFile):
    # Returns a dictionary mapping each protein to its 
    # forward orientation consensus string

    fin = open(origConsensusFile, 'r')
    consensusStr = {}
    for line in fin:
        prot, consensus = line.strip().split()
        consensusStr[prot] = consensus
    return consensusStr

def getOffsets(aliString, forwardCons, A3, pwms):
    # Returns the appropriate offset and rev vectors corresponding to 
    # this alignment (to be used for visualizing 6bp alignment)

    offset, rev = {}, {}
    for p in sorted(aliString.keys()):
        if p not in forwardCons.keys():
            continue
        s = aliString[p]
        f = forwardCons[p]
        k = 0
        for pos in s:
            if pos != '-':
                break
            k += 1
        offset[p] = A3 - k - 2
        hd = 0
        sTrim = s.lstrip('-').rstrip('-')

        # Is it closer to being the forward consensus or the reverse?        
        diff, diff_rc = 0, 0
        f_rc = makeRevCompl(f)
        for i in range(len(f)):
            if sTrim[i] != f[i]:
                diff += 1
            if sTrim[i] != f_rc[i]:
                diff_rc += 1
        if diff_rc < diff:
            f = f_rc
            rev[p] = 1
        else:
            rev[p] = 0
        
        # Just a hack so that creating the offset can't
        # end up outside the range of the PWMs if the alignment
        # is way wrong
        if offset[p] > len(pwms[p]) - 1:
            offset[p] = len(pwms[p]) - 1

        print f[offset[p]:offset[p]+6], rev[p], p
    return offset, rev

def main():

    # Create the output directory
    paramDir = '../STAMP/alignedLogos_NoyesGrps/'
    if not os.path.exists(paramDir):
            os.makedirs(paramDir)

    # Get the data
    if CORE_POS == 'canon9':
        aaPosList = CANON9
    else:
        aaPosList = CORE_POS

    seqs, pwms, core, full, trunc, aaPosList, edges, edges_hmmPos = \
        getHomeoboxData(['n08', 'cisbp', 'chu'], MWID, aaPosList = aaPosList,
                        aposCut = APOS_CUT, edgeCut = EDGE_CUT,
                        maxEdgesPerBase = MAX_EDGES_PER_BASE,
                        N51A_bpos = (MWID-6)/2+2)

    # For mapping/validation
    n08seqs, n08pwms, n08core = seqs['n08'], pwms['n08'], core['n08']
    #cisseqs, cispwms, ciscore = seqs['cisbp'], pwms['cisbp'], core['cisbp']

    #Subset to the same proteins used as input to our procedure
    trSet = 'cisbp'
    pwms1, core1 = pwms[trSet], core[trSet]
    pwms2, core2 = pwms['chu'], core['chu']
    
    #print(len(pwms1), len(pwms2))
    # Remove examples from cis-bp and chu that correspond to 
    # core seqs in the test set
    fin = open('../hd1-explore/0_splitChu2012-trainTest/testProts.txt','r')
    testProts = [x.strip() for x in fin.readlines()]
    fin.close()
    for x in [core1, core2, pwms1, pwms2]:
        allProts = x.keys()
        for prot in allProts:
            if prot in testProts:
                del x[prot]
    #print(len(pwms1), len(pwms2))

    # Combine pwms and core seqs remaining into a single dataset
    pwms, core = {}, {}
    for x in [pwms1, pwms2]:
        for prot in x.keys():
            pwms[prot] = x[prot]
    for x in [core1, core2]:
        for prot in x.keys():
            core[prot] = x[prot]
    #print len(pwms), len(core)

    # Discard proteins with amino acid `X` or '-' in the core sequence
    keep = set([k for k in core.keys() if 'X' not in core[k]])
    [subsetDict(x, keep) for x in [pwms, core]]
    keep = set([k for k in core.keys() if '-' not in core[k]])
    [subsetDict(x, keep) for x in [pwms, core]]
    keep = set([k for k in pwms.keys() if pwms[k].shape[0] >= MWID])
    [subsetDict(x, keep) for x in [pwms, core]]
    
    # Read in the alignment and original consensus string info
    aliString = readSTAMPaliFile(STAMP_ALI_FILE)
    forwardCons = readForwardConsensus(ORIG_CONSENSUS_FILE)

    # Get the starting position and reversal flags for this alignment
    s = aliString['bap']
    f = forwardCons['bap']
    sTrimmed = s.lstrip('-').rstrip('-')
    print s, f, sTrimmed
    # Check for RC alignment
    if f == makeRevCompl(sTrimmed):
        aliString = readSTAMPaliFile(STAMP_ALI_FILE, rc_flag = True)
        s = aliString['bap']
        sTrimmed = s.lstrip('-').rstrip('-')
    assert sTrimmed == f
    A3 = s.index('G')-1
    start, rev = getOffsets(aliString, forwardCons, A3, pwms)

    # Remove anything that wasn't processed by STAMP
    keep = start.keys()
    [subsetDict(x, keep) for x in [pwms, core]]

    print(paramDir)
    # Created the aligned PWMs
    aliPWMS = getAlignedPWMs(getOrientedPWMs(pwms, rev), core, start, 
                             MWID, flipAli = False)

    if False:
        logoDir = paramDir + 'logos_1/'
        if not os.path.exists(logoDir):
            os.makedirs(logoDir)
        makeAllLogos(aliPWMS, core, logoDir)   

    # Compare distance to nearest core/prot in Noyes 08 set
    prot2grp, grp2prot = getNoyesSpecGrps(NOYES_SPEC_GRPS)
    closestCore = getClosestNoyesProt(core, n08core)
    aliScoresClosest = bestAliScoreWithNoyes(closestCore, aliPWMS, n08pwms)

    # Group aligned PWMs that have HD-0 cores with members of Noyes 08 set
    oldPath = paramDir+'logos_1/'
    newPath = paramDir+'logos_grpByNoyes08_hd0/'

    #"""
    if not os.path.exists(newPath):
        os.makedirs(newPath)
    fnames = os.listdir(oldPath)
    fout = open(paramDir+'n08_closestCore.txt', 'w')
    fout.write('prot\tprot.n08\tcore\tcore.n08\tgrp.n08\thDist\tbestAli\n')
    for k in sorted(closestCore.keys()):
        try:
            p, p_n08, dist, c, c_n08, g, aliScore= \
                k, closestCore[k][0], closestCore[k][1], core[k], \
                n08core[closestCore[k][0]], prot2grp[closestCore[k][0]], \
                aliScoresClosest[k]
        except KeyError:
            print "Key not found"
            continue
        fout.write('%s\t%s\t%s\t%s\t%s\t%d\t%f\n' \
                   %(p, p_n08, c, c_n08, g, dist, aliScore))
        if dist == 0:
            oldLab = '_'.join([p,c])
            oldFile = [x for x in fnames if re.match(oldLab,x) is not None][0]
            newFile = '_'.join([g,p_n08])+'_'+oldFile
            os.system('cp %s %s' %(oldPath+oldFile,newPath+newFile))
    fout.close()
    #"""

if __name__ == '__main__':
    main()
