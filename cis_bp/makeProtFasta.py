# Make a fasta file for the protein sequences corresponding to particular 
# motif instances for each TF

PROT_FILE = 'prot_seq.txt'
MOTIF_LIST_FILE = 'motifTable_mostRecent_noMuts.txt'
OUT_FILE = 'prot_mostRecent_noMuts.fa'

# Get the motif IDs of interest
fin = open(MOTIF_LIST_FILE, 'r')
fin.readline()
tfnames = set()
for line in fin:
    tfnames.add(line.strip().split('\t')[0])
fin.close()

# Read in their protein sequences according to CisBP
protSeqs = {}
fin = open(PROT_FILE, 'r')
fin.readline()
for line in fin:
    l = line.strip().split('\t')
    tf = l[2]
    if tf not in tfnames:
        continue
    if ',' not in l[8]:  # Exclude if multiple domains
        if l[9][-1] == '*':
            l[9] = l[9][:-1]
        if tf not in protSeqs.keys():
            protSeqs[tf] = l[9]
        elif tf in protSeqs.keys() and len(l[9]) > len(protSeqs[tf]):
            protSeqs[tf] = l[9]
fin.close()

# Write them to a fasta
fout = open(OUT_FILE, 'w')
for tf in sorted(protSeqs.keys()):
    fout.write('>%s\n' %tf)
    fout.write('%s\n' %protSeqs[tf])
fout.close()