# Make PWM input for the STAMP program to use as a baseline comparison
# for the ability to align the motifs without use of protein information

import numpy as np
from gibbsAlign_GLM import getHomeoboxData 
from getHomeoboxConstructs import subsetDict

BASE = ['A','C','G','T']
B2IND = {x: i for i, x in enumerate(BASE)}
IND2B = {i: x for i, x in enumerate(BASE)}
CORE_POS = 'useStructInfo'    # Uses the precomputed structural alignments
OBS_GRPS = 'grpIDcore'        # Perform group updates based on common "core" AAs in proteins
APOS_CUT = 'cutAApos_1.0_0.05'  # Controls AA contact threshold for binarization
EDGE_CUT = 'edgeCut_1.0_0.05'   # Controls AA contact threshold for binarization
MAX_EDGES_PER_BASE = None       # None means to just ignore this parameter
MWID = 6

OUT_FILE = '../STAMP/inputFile.txt'
OUT_CONSESUS = '../STAMP/stampConsensusCodeLabels.txt'

def getConsensusCode(v):
    # Returns the consensus code consistent with STAMP output
    # given a PWM column, v

    if max(v) >= 0.6:
        return IND2B[np.argmax(v)]
    if v[B2IND['A']]+v[B2IND['C']] >= 0.8:
        return 'M'
    if v[B2IND['A']]+v[B2IND['G']] >= 0.8:
        return 'R'
    if v[B2IND['A']]+v[B2IND['T']] >= 0.8:
        return 'W'
    if v[B2IND['C']]+v[B2IND['G']] >= 0.8:
        return 'S'
    if v[B2IND['C']]+v[B2IND['T']] >= 0.8:
        return 'Y'
    if v[B2IND['G']]+v[B2IND['T']] >= 0.8:
        return 'K'
    return 'N'    

def main():

    # Get the data
    if CORE_POS == 'canon9':
        aaPosList = CANON9
    else:
        aaPosList = CORE_POS

    seqs, pwms, core, full, trunc, aaPosList, edges, edges_hmmPos = \
        getHomeoboxData(['cisbp','chu'], MWID, aaPosList = aaPosList,
                        aposCut = APOS_CUT, edgeCut = EDGE_CUT,
                        maxEdgesPerBase = MAX_EDGES_PER_BASE,
                        N51A_bpos = (MWID-6)/2+2)

    #Subset to the same proteins used as input to our procedure
    trSet = 'cisbp'

    # Subset to the test group and get orientation of key PFM
    pwms1, core1 = pwms[trSet], core[trSet]
    pwms2, core2 = pwms['chu'], core['chu']

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

    # Combine pwms and core seqs remaining into a single dataset
    pwms, core = {}, {}
    for x in [pwms1, pwms2]:
        for prot in x.keys():
            pwms[prot] = x[prot]
    for x in [core1, core2]:
        for prot in x.keys():
            core[prot] = x[prot]

    # Discard proteins with amino acid `X` or '-' in the core sequence
    keep = set([k for k in core.keys() if 'X' not in core[k]])
    [subsetDict(x, keep) for x in [pwms, core]]
    keep = set([k for k in core.keys() if '-' not in core[k]])
    [subsetDict(x, keep) for x in [pwms, core]]
    keep = set([k for k in pwms.keys() if pwms[k].shape[0] >= MWID])
    [subsetDict(x, keep) for x in [pwms, core]]

    fout = open(OUT_FILE, 'w')
    fout_code = open(OUT_CONSESUS, 'w')
    for tf in pwms.keys():
        fout.write('%s  %s  %s\n' %('DE',tf,'Homeobox'))
        m = pwms[tf]
        code = ''
        for i in range(len(m)):
            #bestBase = IND2B[np.argmax(m[i])]
            consensus = getConsensusCode(m[i])
            fout.write('%d  '%i)
            for j in range(len(m[i])):
                fout.write('%f  '%m[i,j])
            fout.write('%s\n'%consensus)
            code += consensus
        fout_code.write('%s\t%s\n' %(tf,code))
        fout.write('XX\n')
    fout.close()
    fout_code.close()

if __name__ == '__main__':
    main()