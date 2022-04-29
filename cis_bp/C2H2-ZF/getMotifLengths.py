# Create a table showing the lengths of each motif

IN_FILE = 'PWM.txt'
OUT_FILE = 'PWM_lengths.txt'

def main():
    fin = open(IN_FILE, 'r')
    fout = open(OUT_FILE, 'w')
    fout.write('Motif_ID\tmotifLen\n')
    line = fin.readline()
    while not line == '':
        l = line.strip().split()
        if len(l) == 0:
            continue
        if l[0] == 'Motif':
            motif = l[1]
        elif l[0] == 'Pos':
            line = fin.readline()
            i = 0
            while line != '\n':
                i += 1
                line = fin.readline()
            fout.write('%s\t%d\n' %(motif, i))
            line = fin.readline()
        line = fin.readline()
    fout.close()


if __name__ == '__main__':
    main()