# Make a fasta file for the protein sequences corresponding to particular 
# motif instances for each TF

PROT_FILE = 'prot_seq_fewZFs.txt'
MOTIF_LIST_FILE = 'motifTable_mostRecent_fewZFs.txt'
OUT_FILE = 'prot_mostRecent_fewZFs.fa'

# Get the motif IDs of interest
fin = open(MOTIF_LIST_FILE, 'r')
fin.readline()
tfnames = set()
for line in fin:
    tfnames.add(line.strip().split('\t')[0])
fin.close()
#print tfnames

# Read in their protein sequences according to CisBP
protSeqs = {}
fin = open(PROT_FILE, 'r')
fin.readline()
for line in fin:
    l = line.strip().split('\t')
    tf = l[0]
    if tf not in tfnames:
        print tf
        continue
    if tf not in protSeqs.keys():
        if l[9][-1] == '*':
            l[9] = l[9][:-1]
        protSeqs[tf] = l[9]
fin.close()

# Write them to a fasta
fout = open(OUT_FILE, 'w')
for tf in sorted(protSeqs.keys()):
    fout.write('>%s\n' %tf)
    fout.write('%s\n' %protSeqs[tf])
fout.close()