# An example showing how to extract interpretable coefficients 
# for per base position amino acid-to-base contact relative energy 
# contributions learned by rCLAMPS


# NOTE:  Set DOMAIN_TYPE in predictionExamples_helpers.py 
#        to 'homeodomain' or 'zf-C2H2' prior to running

from predictionExample import *
OUTFILE = '../coefficientExplore/'+DOMAIN_TYPE+'/coefTable.txt'

def makeCoefTable(model, edges_hmmPos, outfile):
    # Make a table of the model coefficients for exploration
    
    fout = open(outfile, 'w')
    fout.write('bpos\taapos\tbase\taa\tcoef\n')
    for j in range(MWID):
        coefs = model[j].coef_
        intercept = model[j].intercept_
        #print intercept
        aa_pos = edges_hmmPos[j]
        bpos = j+1
        for bInd in range(len(coefs)):
            base = IND2B[bInd]
            #print(len(coefs[bInd]))
            for m in range(len(aa_pos)):
                apos = edges_hmmPos[j][m]
                for aaInd in range(19):
                    aa = IND2A[aaInd]
                    val = coefs[bInd][m*19+aaInd]# + intercept[bInd]
                    fout.write('%d\t%d\t%s\t%s\t%f\n' \
                               %(bpos, apos, base, aa, val))
            fout.write('%d\t%s\t%s\t%s\t%f\n' \
                       %(bpos, 'NA', base, 'NA', intercept[bInd]))

    fout.close()

def main():
    model, aaPosList, edges = createFinalModel(DOMAIN_TYPE)
    
    # Re-index to hmm positions for easier intepretation
    edges_hmmPos = {}
    for bpos in edges.keys():
        edges_hmmPos[bpos] = []
        for apos in edges[bpos]:
            edges_hmmPos[bpos].append(aaPosList[apos])
    makeCoefTable(model, edges_hmmPos, OUTFILE)


if __name__ == '__main__':
    main()

