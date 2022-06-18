# Make PWM input for the STAMP program to use as a baseline comparison
# for the ability to align the motifs without use of protein information

import numpy as np

BASE = ['A','C','G','T']
B2IND = {x: i for i, x in enumerate(BASE)}
IND2B = {i: x for i, x in enumerate(BASE)}
OUT_FILE = '../STAMP/zfC2H2/inputFile.txt'
OUT_CONSESUS = '../STAMP/zfC2H2/stampConsensusCodeLabels.txt'

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

    pwms = readPWMtab('../my_results/zf-C2H2_250_50_seedFFSdiverse6/pwmTab.txt')
    fout = open(OUT_FILE, 'w')
    fout_code = open(OUT_CONSESUS, 'w')
    for tf in pwms.keys():
        fout.write('%s  %s  %s\n' %('DE',tf,'zfC2H2'))
        m = pwms[tf]
        code = ''
        for i in range(len(m)):
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