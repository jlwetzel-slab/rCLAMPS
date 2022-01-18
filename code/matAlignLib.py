# Adjusted functions from Anton's code to do ungapped PWM
# alignment using IC-weighted PCC as the metric

import numpy

def PCC(v1,v2): #compute Pearson Correlation Coefficient for two vectors (positions)
    if len(v1) != 4 or len(v2) != 4:
        sys.stderr.write('WARNING: vectors have wrong length!\n') #must be 4 for nucleotides
        return
    else:
        #get averaged values
        av1, av2 = 0.0, 0.0
        for j in range(4):
            av1 += v1[j]/4
            av2 += v2[j]/4
        #compute PCC
        XY, X2, Y2 = 0.0, 0.0, 0.0
        for j in range(4):
            XY += (v1[j]-av1)*(v2[j]-av2)
            X2 += (v1[j]-av1)**2
            Y2 += (v2[j]-av2)**2
        if X2*Y2 != 0:
            pc = XY / (X2*Y2)**0.5
        else: #correction for the vectors like (0.25,0.25,0.25,0.25)
            pc = 0.0
    return pc

def information_content(column): #compute information content of the column
    E = 0.0 #entropy
    for j in range(4):
        if column[j] != 0:
            try:
                E -= column[j] * numpy.log2(column[j])
            except RuntimeWarning:
                print column[j]
    ic = 2 - E
    return ic

def norm_matrix(matrix):
    #normalize each column to 1
    columns = len(matrix) #get number of columns
    rows = len(matrix[0]) #get number of rows
    matrix_n = numpy.zeros((columns,rows), float)
    for i in range(columns):
        #read column
        summa = 0
        for j in range(rows):
            summa += matrix[i][j]
        #normalize column
        for j in range(rows):
            matrix_n[i][j] = float(matrix[i][j]) / summa
    return matrix_n

def matrix_align(mp, me, metric, normScore = False, oneSided = False, minWidthM2 = 1):
    #align two matrices probing all possible shifts (without Smith-Waterman algorithm)
    #score takes into account the IC content of both columns symmetrically
    L1 = len(mp)
    L2 = len(me)
    #normalize matrices)
    mp = norm_matrix(mp)
    me = norm_matrix(me)

    best_score, best_shift = 0.0, 0
    best_ic_only = 0.0
    #try all possible shifts
    for shift in range(minWidthM2-L2, L1-minWidthM2+1):
        i,j = 0,0 #starting position
        if shift < 0: #beginning of sequence 2 is cut
            j = -shift
        elif shift > 0: #beginning of the sequence 1 is cut
            i = shift
        #get the diagonal score
        score, ic_only = 0.0, 0.0
        while i<L1 and j<L2:
            m = eval(metric)(mp[i],me[j])
            if oneSided:
                ic_norm = information_content(me[j])/2
            else:
                ic_norm = (information_content(mp[i])+information_content(me[j]))/4
            score += m*ic_norm
            ic_only += ic_norm
            i+=1; j+=1
        if score > best_score:
            best_score = score
            best_shift = shift
            best_ic_only = ic_only
        #if ic_only >= best_ic_only:
        #    best_ic_only = ic_only

    if normScore:
        best_score /= best_ic_only
    return (best_score, best_shift)

def matrix_compl(matrix):
    #get reversed complement of the matrix
    m_comp = numpy.zeros((len(matrix),4), float)
    for i in range(len(matrix)):
        for j in range(4):
            m_comp[len(matrix)-i-1][3-j] = matrix[i][j]
    return(m_comp)

def comp_matrices(mp, me, metric = 'PCC', normScore = False,
                  useMaxInstead = False, getAvgPCC = False,
                  oneSided = False, minWidthM2 = 1):
    #compare predicted PWM against experimental PWM
    score, shift = matrix_align(mp, me, 'PCC', normScore = normScore,
                                oneSided = oneSided,
                                minWidthM2 = minWidthM2)
    score2, shift2 = matrix_align(mp, matrix_compl(me), 'PCC',
                                  normScore = normScore, 
                                  oneSided = oneSided, 
                                  minWidthM2 = minWidthM2)
    if score2 > score:
        score = score2
        shift = shift2
        rev = True #reversed flag
    else:
        rev = False

    #### For testing purposes #####
    if useMaxInstead:
        x,_ = matrix_align(mp, mp, 'PCC', normScore = False)
        y,_ = matrix_align(me, me, 'PCC', normScore = False)
        score /= max(x, y)

    # Return the average (uncorrected) PCC per column as well
    #if getUncorrectedScore:
    #    avgPCC = 0.0
    #    if rev:
    #        me = matrix_compl(me)
    #    i,j = 0,0 #starting position
    #    if shift < 0: #beginning of sequence 2 is cut
    #        j = -shift
    #    elif shift > 0: #beginning of the sequence 1 is cut
    #        i = shift


    #    return (score shift, rev, avgPCC)

    return (score,shift,rev)
