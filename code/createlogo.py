import pickle
import numpy as np
from matAlignLib import *
from pwm import makeNucMatFile, makeLogo
from copy import deepcopy
import os

MWID = 8

def getOrientedPWMs(pwms, rev, reorient = 1):
    """
    Algorithm function.
    Returns PWMs according to rev
    """

    # Orient pwms according to rev
    pwms_o = {}
    for k in pwms.keys():
        if rev[k] == reorient:
            pwms_o[k] = matrix_compl(pwms[k])
        else:
            pwms_o[k] = deepcopy(pwms[k])
    return pwms_o

def makeAllLogos(pwms, aSeqs, logoDir, keysToUse = None):
    """
    Output function
    Place logos for every pwm and ouput to logoDir
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

def getAlignedPWMs(pwms, aSeqs, start, mWid, flipAli = False):
    """
    Output function.
    Returns a new set of PWMs, truncated on each side to the
    aligned region
    """

    npwms = {}
    for p in pwms.keys():
        pwm = pwms[p]
        npwm = np.zeros((mWid,4), dtype = 'float')
        s = start[p]
        for i in range(mWid):
            npwm[i,:] = pwm[i+s,:]
        if flipAli:
            npwm = matrix_compl(npwm)
        npwms[p] = npwm
    return npwms

def main():
    mWid = MWID
    print("Creating aligned logos ...")

    resfile = open("res.obj",'rb')
    res = pickle.load(resfile)
    resfile.close()

    mdfile = open("md.obj",'rb')
    mainOutDir = pickle.load(mdfile)
    mdfile.close()

    pwmsfile = open("pwms.obj",'rb')
    pwms = pickle.load(pwmsfile)
    pwmsfile.close()

    corefile = open("core.obj",'rb')
    core = pickle.load(corefile)
    corefile.close()

    ll = [x['ll'] for x in res]
    S = [x['S'] for x in res]
    O = [x['O'] for x in res]
    reorient = [x['reorient'] for x in res]
    opt = np.argmax(ll)


    for i, o in enumerate([opt]):
        logoDir = mainOutDir + 'logos_%d/' %(i+1)
        if not os.path.exists(logoDir):
            os.makedirs(logoDir)
        flipAli = False
        if reorient[o]:
            flipAli = True
        aliPWMs = getAlignedPWMs(getOrientedPWMs(pwms, O[o]),
                                    core, S[o], mWid, flipAli = flipAli)
        makeAllLogos(aliPWMs, core, logoDir)


if __name__ == '__main__':
    main()

