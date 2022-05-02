import os
import numpy as np
from getHomeoboxConstructs import readFromFasta, writeToFasta, subsetDict
# parse a txt file downloaded from cis-bp

CISBP_SEQ_fpath_txt = "../cis_bp/prot_seq.txt"
CISBP_PWM_fpath = "../cis_bp/PWM.txt"

# fa file form:
# >TF_ID
# protein sequence
def convertToFasta(fpath_txt):
    fin = open(fpath_txt)
    line = fin.readline()
    fpath_fa = os.path.splitext(fpath_txt)[0]+'.fa'
    file = open(fpath_fa, 'w+')
    file.truncate(0)
    while line != '':
        line = fin.readline()
        lineArr = line.split('\t')
        if len(lineArr) < 10:
            continue
        if len(lineArr[8].split(",")) > 1:
            continue
        file.write(">"+lineArr[2]+"\n")
        file.write(lineArr[8]+"\n")
    file.close()
    return fpath_fa

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

def getPWM_barrera(fpath, motifs, tfnFull, verbose = False):
    # Extracts motifs from the CIS-BP PWM.txt file subset to 
    # only those whose TF_Name is in tfs (set) and whose 
    # Motif_ID is in motifs (set)
    fin = open(fpath)
    pwms = {}
    line = fin.readline()

    while line != "":
        line = fin.readline()
        lineArr = line.split("\t")
        #if verbose:
        #    print lineArr
        if lineArr[0] == "TF Name":
            tf = lineArr[1].rstrip()
        if lineArr[0] == "Pos":
            pwm = []
        if lineArr[0] == "Motif":
            motif = lineArr[1].rstrip()
        if len(lineArr) == 5 and lineArr[0] != "Pos":
            lineArr = np.array(lineArr)
            onevec_list = lineArr.astype(np.float)[1:5].tolist()
            pwm.append(onevec_list)
        if lineArr[0] == '\n' and motif in motifs:
            pwms[tfnFull[motifs.index(motif)]] = np.array(pwm)
    return pwms




def main():
    seqs = {}   # The full length protein construct
    pwms = {}

    dset = 'cisbp'
    CISBP_FASTA = convertToFasta("../cis_bp/prot_seq.txt")
    CISBP_PWM_SRC = "../cis_bp/PWM.txt"
    seqs[dset] = readFromFasta(CISBP_FASTA)
    pwms[dset] = getPWM(CISBP_PWM_SRC)
            # [FIXME]: not sure about this subset

    subsetDict(pwms[dset], set(seqs[dset].keys()))
    stem = '/'.join(CISBP_FASTA.split('/')[:-1]) + \
            '/homeodomains_cisbp'

    #convertToFasta(CISBP_SEQ_fpath_txt)
    #getPWM(CISBP_PWM_fpath)


if __name__ == '__main__':
    main()

