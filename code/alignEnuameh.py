# A script for determining the correct orientation of the
# full-length Enuameh PWMs by finding the optimal alignment of 
# each to its single-finger concatenated version going in reverse-finger order

import numpy as np
from matAlignLib import comp_matrices
from gibbsAlign_GLM import getPrecomputedInputs_zfC2H2

def getPWMsRevFingOrder(fname):
    fin = open(fname, 'r')
    line = fin.readline()
    line = fin.readline()
    currentProt = ''
    pwms = {}
    i = 0
    while line != '':
        #print i
        l = line.strip().split('\t')
        lastProt = currentProt
        currentProt = l[0]
        if currentProt != lastProt:
            if i != 0:
                pwms[lastProt] = pwm
            pwm = []
        pwm.append([float(x) for x in l[4:]])
        line = fin.readline()
        i += 1
    pwms[lastProt] = pwm
    fin.close()

    return pwms

def main():
    print "Hello"
    pwms_full, _, _, _, _ = \
        getPrecomputedInputs_zfC2H2(rescalePWMs = False, ffsOnly = True, includeB1H = False)
    pwms_revFingOrder = getPWMsRevFingOrder('../flyFactorSurvey/enuameh/enuameh_perFinger_PWMs_reverseFingerOrder.txt')
    fout = open('../flyFactorSurvey/enuameh/enuameh_startPosInfo.txt', 'w')
    fout.write('prot\tstart\trev\n')
    for p in sorted(pwms_full.keys()):
        #print len(pwms_full[p]), len(pwms_revFingOrder[p])
        if len(pwms_revFingOrder[p]) > len(pwms_full[p]):
            continue
        pWid = len(pwms_full[p])
        
        # Align the annotated finger PWM to the 
        score, shift, ori = comp_matrices(pwms_revFingOrder[p], pwms_full[p])
        #score, shift, ori = comp_matrices(pwms_revFingOrder[p], pwms_full[p], 
        #                                  minWidthM2 = len(pwms_revFingOrder[p]), 
        #                                  oneSided=True)
        if p == 'CG10267':
            print pwms_full[p]
            print pwms_revFingOrder[p]
            print p, score, -shift, int(ori)
        fout.write('%s\t%d\t%d\n' %(p,-shift,int(ori)))
    fout.close()


if __name__ == '__main__':
    main()