# Make a DBD sequence FASTA file for the Chu et al. mutants
# by substituting the correct positions within the DBD for 
# the Engrailed protein (en) as given in protein sequences from 
# the Cis_BP database

MATCH_TAB = '../cis_bp/homeodomains_cisbp_hasPWM.matchTab.txt'
HMM_POS = [43,46,47,50,54]#[2,3,4,5,47,50,51,54,55]#
PWM_TAB = 'MEME_allMotifs_rw.txt'
SUBSET_TAB = 'allMotifs_info.txt'  # Just outputs a simple file with all motif listings

def main():

    # Read in the en protein
    fin = open(MATCH_TAB, 'r')
    for line in fin:
        if line[:3] == 'en ':
            en = ''.join(line.strip().split()[1:])
            break
    fin.close()

    # Read in the set of TF names from the PWM table
    fin = open(PWM_TAB, 'r')
    tfNames = []
    motifs = []
    for line in fin:
        if line[:7] == 'TF Name':
            tfNames.append(line[8:].rstrip())
        if line[:5] == 'Motif':
            motifs.append(line[6:].rstrip())
    fin.close()

    # Create the fasta file
    subPos = [x-2 for x in HMM_POS]
    fout = open('domains.fasta', 'w')
    for tfn in tfNames:
        fout.write('>%s\n' %tfn)
        newprot = ''
        subs = tfn.split('_')[-1]
        #print subs
        i = 0
        for k, aa in enumerate(en):
            if i < len(subPos) and k == subPos[i]:
                newprot += subs[i]
                i += 1
            else:
                newprot += aa
        fout.write(newprot+'\n')
    fout.close() 

    # Create the motif listings table
    fout = open(SUBSET_TAB, 'w')
    for i in range(len(tfNames)):
        fout.write('%s\t%s\n' %(tfNames[i], motifs[i]))
    fout.close()

if __name__ == '__main__':
    main()


