# A script for making predictions based on parameter
# estimates output by gibbAlign_GLM.py [predict fly proteins using the model built using mouse proteins]

from gibbsAlign_GLM import getPrecomputedInputs_zfC2H2, readSeedAlignment
from gibbsAlign_GLM import makeAllLogos, assignObsGrps
from gibbsAlign_GLM import getAlignedPWMs_multiDomain, getOrientedPWMs
from getHomeoboxConstructs import subsetDict
import numpy as np
import pickle, os

BASE = ['A','C','G','T']
REV_COMPL = {'A':'T','C':'G','G':'C','T':'A'}
AMINO = ['A','C','D','E','F','G','H','I','K','L',
         'M','N','P','Q','R','S','T','V','W','Y']
B2IND = {x: i for i, x in enumerate(BASE)}
A2IND = {x: i for i, x in enumerate(AMINO)}
IND2B = {i: x for i, x in enumerate(BASE)}
IND2A = {i: x for i, x in enumerate(AMINO)}

# Input files for PWMs and protein info
PROT_SEQ_FILE = '../precomputedInputs/zf-C2H2/prot_seq_fewZFs_hmmerOut_clusteredOnly_removeTooShort.txt'  # Input protein domain file subsetted to relvant amino acid contacting positions
PROT_SEQ_FILE_FFS = '../flyFactorSurvey/enuameh/enuameh_perFinger_processedProtInfo.txt'
PWM_INPUT_TABLE = '../precomputedInputs/zf-C2H2/pwmTab_fewZFs_clusteredOnly_removeTooShort.txt'   # A table of PWMs corresponding to prots in PROT_SEQ_FILE
PWM_INPUT_FILE_FFS = '../flyFactorSurvey/enuameh/flyfactor_dataset_A.txt'
SEED_FILE = '../flyFactorSurvey/enuameh/enuameh_startPosInfo.txt'

OBS_GRPS = 'grpIDcore'
MWID = 4
RIGHT_OLAP = 1
ANCHOR_B1H = False

TAB_DIR = '../my_results/zf-C2H2_250_50_seedFFSdiverse6/'
MODEL_FILE = TAB_DIR+'result.pickle'
OUT_DIR = TAB_DIR+'plots/'

def main():

    with open(MODEL_FILE) as f:
        res = pickle.load(f)

    score = [x['ll'] for x in res]
    reorient = [x['reorient'] for x in res]
    start = [x['start'] for x in res]
    rev = [x['rev'] for x in res]
    opt = np.argmax(score)
    print(opt)


    pwms, core, edges, edges_hmmPos, aaPosList = \
        getPrecomputedInputs_zfC2H2(rescalePWMs=False,ffsOnly=False,includeB1H=False)
    start = start[opt]
    rev = rev[opt]
    
    # Remove examples where PWMs that are too short for the number of domains
    nDoms = {}
    for p in core.keys():
        nDoms[p] = len(core[p])/len(aaPosList)
        if len(pwms[p]) < (MWID-RIGHT_OLAP)*nDoms[p]+RIGHT_OLAP or nDoms[p] < 2:
            del nDoms[p]
            del core[p]
            del pwms[p]
        
    # Remove examples where the known/stated fixed starting position
    # would make the pwm too short for the number of arrays annotated as binding
    knownStarts_ffs = readSeedAlignment(SEED_FILE, include = pwms.keys())
    for p in knownStarts_ffs.keys():
        if len(pwms[p]) < knownStarts_ffs[p]['start']+(MWID-RIGHT_OLAP)*nDoms[p]+RIGHT_OLAP:
            #print p, core[p], nDoms[p], len(pwms[p])
            del nDoms[p]
            del core[p]
            del pwms[p]
            del knownStarts_ffs[p]

    # Assign to observation groups with identical core sequences
    obsGrps = assignObsGrps(core, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots  

    # Evaluate the fit of the aligned pwms to the inferred model
    flipAli = False
    if reorient[opt]:
        flipAli = True

    # Get the PWM alignments inferred for the non-B1H, multi-finger ZF arrays
    fout = open(TAB_DIR+'registrationInfo.txt', 'w')
    fout.write('prot\tstart\trev\n')
    for p in sorted(start.keys()):
        if p[:4] != 'B1H.':
            #print p,start[p], rev[p]
            fout.write('%s\t%d\t%d\n' %(p,start[p],rev[p]))
    fout.close()

    #"""
    aliPWMS = getAlignedPWMs_multiDomain(getOrientedPWMs(pwms, rev), core, start,
                                         nDoms, flipAli = flipAli)
    #for k in aliPWMS.keys():
    #    print k, len(aliPWMS[k])
    #    print aliPWMS[k]
    logoDir = OUT_DIR + '0_logos_aligned/'
    print "Creating aligned logos in %s" %logoDir
    if not os.path.exists(logoDir):
        os.makedirs(logoDir)
    makeAllLogos(aliPWMS, core, logoDir)
    #"""


if __name__ == '__main__':
    main()
