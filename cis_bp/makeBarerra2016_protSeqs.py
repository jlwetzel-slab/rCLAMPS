# Create the protein sequences corresponding to the Barerra 2016 mutant 
# homeodoman PBMs based on the original protein from the cisBP file

ORIG_PROT_TAB = 'prot_seq.txt'
BARERRA_TAB = 'motifTable_Barrera2016_mutsOnly.txt'
FASTA_OUT = 'prot_seq_Barrera2016_mutsOnly.fasta'
TAB_OUT = 'mutInfoTable_Barrera2016_mutsOnly.txt'

def getOrigProtSeq(protFile_cisBP, tfname, fr, pos):
    # Returns the full length protein sequence from the 
    # cisBP protein info file corresponding to the tf tfname,
    # checking that this isoform has the correct protien at the 
    # substituiton position

    fin = open(protFile_cisBP, 'r')
    fin.readline()
    #print tfname
    full = ''
    for line in fin:
        l = line.strip().split('\t')
        #print l
        if l[2] == tfname:
            newfull = l[9]
            if len(full) < len(newfull) and newfull[pos-1] == fr:
                full = newfull
                ensg, ensp = l[3], l[4]
    fin.close()
    return full, ensg, ensp


def main():
    
    # Read the table to get the original protein names and the
    # mutation information
    fin = open(BARERRA_TAB, 'r')
    fin.readline()
    foutTab = open(TAB_OUT,'w')
    foutTab.write('TF Name\tGene ID\tProtein ID\tmutPos\tfrom\tto\tmotif\n')
    fastaOut = open(FASTA_OUT,'w')
    for line in fin:
        l = line.strip().split('\t')
        tf, mInfo, motifNew = l[0], l[13], l[5]
        mutInfo = mInfo.split('_')[1]
        fr, pos, to = mutInfo[0], int(mutInfo[1:-1]), mutInfo[-1]
        full, ensg, ensp = getOrigProtSeq(ORIG_PROT_TAB,tf, fr, pos)
        #print tf, full
        if full == '':
            continue
        assert fr == full[pos-1]
        foutTab.write('%s\t%s\t%s\t%d\t%s\t%s\t%s\n' %(mInfo, ensg, ensp, pos, fr, to, motifNew))
        newprot = full[:pos-1]+to+full[pos:]
        assert len(newprot) == len(full)
        fastaOut.write('>%s\n%s\n' %(mInfo, newprot))
    foutTab.close()
    fastaOut.close()
    fin.close()

if __name__ == '__main__':
    main()

