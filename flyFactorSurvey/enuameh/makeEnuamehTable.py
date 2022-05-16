# Create a table from the Enuameh supplemental data including all
# of the SOLEXA motifs annotated with finger assignments

IN_FILE = 'Supplemental_Dataset_1.txt'
OUT_FILE = 'enuameh_perFinger.txt'

def main():
    fin = open(IN_FILE, 'r')
    fout = open(OUT_FILE, 'w')
    fout.write('prot\tzfNum\thelix\tcore\tmotif\tmotifPos\tstrand\tA\tC\tG\tT\n')
    line = fin.readline()
    while line != '':
        assert line[0] == '>'
        l = line.strip().split('\t')
        helix = l[0][1:]
        print helix
        core = helix[0]+helix[2]+helix[3]+helix[5]        
        zfNum = l[3][1:]
        prot = l[4].split('_')[0]
        motif = l[4].split('.')[0]
        strand = l[5].split('=')[1]
        motifPos = int(l[6].split('=')[1])
        if 'SOLEXA' not in motif:
            [fin.readline() for k in range(4)]
            line = fin.readline()
            continue
        pwm = [['','','',''] for k in range(4)]
        print [prot,zfNum,str(motifPos),strand]
        for i in range(4):
            line = fin.readline()
            l = line.strip().split('\t')
            for j in range(len(l)):
                pwm[j][i] = l[j]
        for i in range(len(pwm)):
            fout.write('\t'.join([prot,zfNum,helix,core,motif,str(motifPos+i),strand]+pwm[i])+'\n')
        line = fin.readline()
    fin.close()
    fout.close()

if __name__ == '__main__':
    main()