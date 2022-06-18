# Parse the output of STAMP to see if any offset exists that
# would result in a good alignment of the C2H2-ZF PWMs for 
# which we know the starting positions (from Enuameh 2013)

import numpy as np
import os, re
#from gibbsAlign_naiveBayes import getAlignedPWMs, getOrientedPWMs, reverseBaseOrient
#from getHomeoboxConstructs import subsetDict

STAMP_ALI_FILE = '../STAMP/zfC2H2/Stamp_Results_zfC2H2_enuamehOnly.txt'
ORIG_CONSENSUS_FILE = '../STAMP/zfC2H2/stampConsensusCodeLabels.txt'
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

def readActualStarts(startFile):
    # Read in the actual starting positions for each 
    # of the enuameh PWMs
    starts = {}
    fin = open(startFile)
    fin.readline()
    for line in fin:
        l = line.strip().split()
        starts[l[0]] = {'start': int(l[1]), 'rev': int(l[2])}
    fin.close()
    return starts

def main():
    # Read in the alignment and original consensus string info
    aliString = readSTAMPaliFile(STAMP_ALI_FILE)
    aliString_rev = {k: makeRevCompl(aliString[k]) for k in aliString.keys()}
    forwardCons = readForwardConsensus(ORIG_CONSENSUS_FILE)
    forwardCons_rev = {k: makeRevCompl(forwardCons[k]) for k in forwardCons.keys()}

    starts = readActualStarts('../flyFactorSurvey/enuameh/enuameh_startPosInfo.txt')
    for k in aliString.keys():
        if k not in starts.keys():
            del aliString[k]
            del aliString_rev[k]
        if k in forwardCons.keys():
            del forwardCons[k]

    # For each possible TF, assume it is aligned correctly and
    # count how many others are aligned correctly under that assumption
    maxCorrect = 1
    for p in aliString.keys():
        if starts[p]['rev'] == 1:
            a = aliString_rev
            f = forwardCons_rev
        else:
            a = aliString
            f = forwardCons
        offset = 0
        for x in a[p]:
            if x != '-':
                break
            else:
                offset += 1
        thisStart = offset+starts[p]['start']
        nCorrect = 1
        for p2 in aliString.keys():
            if p == p2:
                continue
            offset2 = 0
            for x in a[p2]:
                if x != '-':
                    break
                else:
                    offset2 += 1
            thisStart2 = offset2+starts[p2]['start']
            if thisStart == thisStart2:
                nCorrect += 1
        if nCorrect > maxCorrect:
            maxCorrect = nCorrect

    fracCorrect = nCorrect/float(len(aliString.keys()))*100
    print "At most %f%% of C2H2-ZF mappings are correct in the STAMP alignment." %fracCorrect


if __name__ == '__main__':
    main()
