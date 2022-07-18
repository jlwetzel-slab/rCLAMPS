# Uses runhmmer.py functions to gather a list of domains
# for a given set of database files.

from runhmmer import runhmmer232, getdescs, idgen, sorttable
import os


def checkProperties(seq, aliseq):
    # Checks for characteristic properties of a zinc
    # finger sequence with respect to alignment to the
    # PFAM HMM, returning true if the properties
    # hold, false otherwise.  seq should be the model
    # string output by hmmer2.3.2 and aliseq should
    # be the portion of sequence aligned to the model

    aliseq = aliseq.upper()

    # Position of histidines/histidine-cysteine
    if not (seq[-1] == 'H' and (aliseq[-1] == 'H' or aliseq[-1] == 'C')):
        return False, 'noH2'
    H2pos = len(seq) - 1
    H1pos = seq[:H2pos].rfind('H')
    if not aliseq[H1pos] == 'H':
        return False, 'noH1'
    if not H2pos-H1pos >= 3 and H2pos-H1pos <= 6:
        return False, 'HpairDist'

    # Position of cysteines
    C1pos = seq.find('C')
    if not (C1pos == 2 and aliseq[C1pos] == 'C'):
        return False, 'noC1'
    C2pos = seq[3:H1pos].find('C') + 3
    if not aliseq[C2pos] == 'C':
        return False, 'noC2'
    if not C2pos-C1pos >= 3 and C2pos-C1pos <= 6:
        return False, 'CpairDist'

    # Require exactly 12 residues between C2pos and H1pos
    if '.' in seq[H1pos-7:H1pos]:
        return False, 'helixIns'
    if '-' in aliseq[H1pos-7:H1pos]:
        return False, 'helixGap'

    return True, None

def getLinkerLens(table, protCol, startCol, endCol, maxLinker = 8):
    # Helper function for reformatTable.  Computes the
    # lengths of the left and right linkers for each
    # domain and assigns domains to clusters.
    # Returns a tuple of lists for leftLinkers
    # rightLinkers, and cluster numbers

    leftLinkers = []
    clusts = []

    currProt = None
    currEnd = table[0][endCol]
    for row in table:
        nextProt = row[protCol]
        nextStart = row[startCol]

        # Decide if we are in a new cluster or not
        # and compute left linker lengths
        if nextProt == currProt:
            leftLinkerLen = (nextStart - currEnd - 1)
            if leftLinkerLen > maxLinker:
                clustNum += 1
                leftLinkerLen = 'NA'
        else:
            clustNum = 0
            leftLinkerLen = 'NA'

        # Append and update variables
        leftLinkers.append(leftLinkerLen)
        clusts.append(clustNum)
        currProt = nextProt
        currEnd = row[endCol]

    # Right linker is left linker of next domain
    # unless that domain is from another cluster,
    # in which case right linker is NA
    rightLinkers = []
    for i in range(len(leftLinkers) - 1):
        currProt = table[i][1]
        currClust = clusts[i]
        nextProt = table[i+1][1]
        nextClust = clusts[i+1]
        if currProt == nextProt and currClust == nextClust:
            rightLinkers.append(leftLinkers[i+1])
        else:
            rightLinkers.append('NA')
    rightLinkers.append('NA')

    return leftLinkers, rightLinkers, clusts

def makeTable(fname, db, checkQual = False):
    # Takes transformed hmmer232 output as input.
    # Retuns two lists:
    # The first is a list of tuples where each element of the
    # tuple corresponds to a column of a table.  The second
    # is a list of names for those columns.

    fin = open(fname)
    fin.readline()

    table = []
    discarded = {x: 0 for x in ['noH1', 'noH2', 'noC1', 'noC2',
        'HpairDist', 'CpairDist', 'helixIns', 'helixGap']}
    
    names = ['prot','targStart','targEnd','helixStart','helixSeq','coreSeq', \
        'hmm','modelSeq','aliSeq','bitscore','e_val']
    for i, line in enumerate(fin):

        l = line.strip().split()
        tfname = l[0]
        hmm = l[1]
        e_val = eval(l[2])
        bitscore = eval(l[3])
        targStart = eval(l[4])
        targEnd = eval(l[5])
        modelSeq = l[6]
        aliSeq = l[7]

        # Check that the alignment looks reasonable
        zfGood, reason = checkProperties(modelSeq, aliSeq)
        if not zfGood:
            discarded[reason] += 1
            continue

        # Compute the helix and core seqs based on starting H
        helixStart = modelSeq.find('H') - 7
        helixSeq = aliSeq[helixStart:helixStart+7]
        coreHelixSeq = aliSeq[helixStart] + aliSeq[helixStart+2] + \
            aliSeq[helixStart+3] + aliSeq[helixStart+6]
        helixStart = helixStart + targStart

        row = [tfname, targStart, targEnd, helixStart, helixSeq, \
            coreHelixSeq, hmm, modelSeq, aliSeq, bitscore, e_val]
        table.append(row)
    fin.close()

    # Sort on gene, pname, pname2, targStart (or the like)
    table = sorttable(table, [0,1])
    return table, names, discarded


def reformatOutput(fname, db = 'uniprot'):
    # Reformat the output that Shilpa's runhmmer232
    # code outputs

    table, names, discarded = makeTable(fname, db = db)
    #print(table)

    # Print basic summary
    print "Found %d valid domains." %len(table)
    print "Domains discarded:"
    for k in discarded.keys():
        print "-- %s:\t%d" %(k, discarded[k])

    #print names
    #print table[0]
    # Get "linker" lengths
    protCol, startCol, endCol = 0, 1, 2
    leftLinkers, rightLinkers, clusts = \
            getLinkerLens(table, protCol, startCol,endCol)

    # Put linker lenghts and cluster numbers into table
    for i, row in enumerate(table):
        table[i] = row + [clusts[i], leftLinkers[i],rightLinkers[i]]
    names = names + ['clustNum', 'leftLinkerLen', 'rightLinkerLen']

    # Write the table to a file
    newFile = './tmp' + idgen() + '.txt'
    fout = open(newFile, 'w')
    fout.write('\t'.join(names) + '\n')
    for row in table:
        fout.write('\t'.join([str(x) for x in row]) + '\n')
    fout.close()
    os.system('mv %s %s' %(newFile, fname))

def runHmmer232OnFasta(inFile, outFile, hmmLab, db = 'uniprot'):
    
    # Runs hmmer232 on a uniprot fasta file.
    descs = getdescs(inFile)
    hmmFile = hmmLab+'.232.hmm'
    runhmmer232(outFile, hmmFile, inFile, hmmLab, descs)
    reformatOutput(outFile)

def make4posTxt(hmmerOutFile, txtFile, colNum):
    # Reads the hmmerOutFile colNum, subsets to the unique
    # entries in this column them to a plain text file
    # one line per entry.

    fin = open(hmmerOutFile, 'r')
    fin.readline() # skip the header
    seqs = set([x.strip().split('\t')[colNum] for x in fin])
    fin.close()
    fout = open(txtFile, 'w')
    for i, seq in enumerate(sorted(list(seqs))):
        fout.write('%s\n' %(seq))
    fout.close()

def main():

    inFile = 'prot_mostRecent_fewZFs.fa'
    outFile = 'prot_seq_fewZFs_hmmerOut.txt'
    runHmmer232OnFasta(inFile, outFile, 'zf-C2H2')
    #fastaFile = '../hmmerOut/%s/%s/%s-%s_4pos.fasta' %(HMM_LAB,DB,DB_SUB,ORG)
    #txtFile = '../hmmerOut/%s/%s/%s-%s_4pos.txt' %(HMM_LAB,DB,DB_SUB,ORG)
    #makeFakeDom4posFasta(outFile, fastaFile, 9)
    #make4posTxt(outFile, txtFile, 9)



if __name__ == '__main__':
    main()

