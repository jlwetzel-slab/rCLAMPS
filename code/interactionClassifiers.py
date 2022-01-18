# Make classifiers with interaction terms based on the PWMs that were 
# alligned by the training procedure

from gibbsAlign_GLM import *
from hd1_transfer_predictions import getTrainPairsAndInfo, getTestProtsAndPWMs
from sklearn.ensemble import RandomForestClassifier
from makeAlignedPWMs import makePCCtable
from sklearn.model_selection import GroupKFold, GridSearchCV, ParameterGrid
from sklearn.metrics import log_loss


#############
## Simply set these params the same as for the Gibbs sampler so 
## the code knows where to find the PWM alignment info
#############

# These are parameters for adjusting the model topology
# (i.e., the binarized contact matrix)
CORE_POS = 'useStructInfo'    # Uses the precomputed structural alignments
OBS_GRPS = 'grpIDcore'        # Perform group updates based on common "core" AAs in proteins
APOS_CUT = 'cutAApos_1.0_0.05'  # Controls AA contact threshold for binarization
EDGE_CUT = 'edgeCut_1.0_0.05'   # Controls AA contact threshold for binarization
MAX_EDGES_PER_BASE = None       # None means to just ignore this parameter

TR_SET = 'cisbp'#'b08' #  # Just controls where code looks for the PCMs

MWID = 6  # Number of base positions in the PDSIM
RAND_SEED = 382738375 #78374223 # Random seed for reproducibilty
MAXITER = 15
N_CHAINS = 100
INIT_ORACLE = False
SAMPLE = 100
MAKE_LOGOS = True

################

HMM_OFFSET = 2
MIP_CORE_POS = [3,6,19,47,50,54,55]
MIP_CORE_POS = [x - HMM_OFFSET for x in MIP_CORE_POS]
FIT_MODEL = True
USE_STRUCT_EDGES = False
DEFAULT_PARAMS = True
USE_MSE_CV = False

"""
def createRandomForestModel(trainX, trainY, trainW):
    #
    #Function for performing eight weighted multinomial logistic regressions for eight base positions
    #:param trainX: a dictionary of {base position (0-7): X matrix}
    #:param trainY: a dictionary of {base position (0-7): y vector}
    #:param trainW: a dictionary of {base position (0-7): weights, same length as y vector}
    #:return:
    #model: a dictionary of {base position (0-7): GLM model}
    # Standarize features [Question]: I'm not sure that whether I should standardize X since X is very sparse
    #scaler = StandardScaler()
    #X_std = scaler.fit_transform(trainX)
    model = {}
    for j in range(MWID):
        #clf = LogisticRegression(fit_intercept=True, random_state=0, multi_class='multinomial', solver='newton-cg')
        clf = DecisionTreeClassifier()
        #clf = DecisionTreeClassifier(max_depth = 15)
        model[j] = clf.fit(trainX[j], trainY[j], sample_weight=trainW[j])
    return model
"""

def grid_cv(X_in, y_in, w_in, cv, param_grid, 
            use_weighting = False, groups = None):
    out_results = dict()

    for i, k in enumerate(param_grid):
        clf = RandomForestClassifier(n_estimators=100,criterion="entropy",
                                     warm_start=False,n_jobs=-1,**k)
        for train_ndx, test_ndx in cv.split(X=X_in, y=y_in, groups = groups):
            X_train = X_in[train_ndx, :]
            y_train = y_in[train_ndx]
            w_train = w_in[train_ndx]
            y_test = y_in[test_ndx]

            clf.fit(X=X_train,y=y_train,sample_weight=w_train)

            y_hat = clf.predict_proba(X=X_in[test_ndx, :])
            if use_weighting:
                w_test = w_in[test_ndx]
                w_i_sum = w_test.sum()
                score = w_i_sum/w_in.sum()*log_loss(y_true=y_test, y_pred=y_hat, sample_weight=w_test)
            else:
                score = log_loss(y_true=y_test,y_pred=y_hat)

            results = out_results.get(i, [])
            results.append(score)
            out_results.update({i: results})

    for k, v in out_results.items():
        if use_weighting:
            mean_score = sum(v)
        else:
            mean_score = np.mean(v)
        out_results.update({k: mean_score})

    best_score = min(out_results.values())
    best_param = min(out_results, key=out_results.get)
    
    return best_score, best_param

def grid_cv_mse(X_in, y_in, w_in, cv, param_grid, groups = None):
    out_results = dict()

    for i, k in enumerate(param_grid):
        clf = RandomForestClassifier(n_estimators=100,criterion="entropy",
                                     warm_start=False,n_jobs=-1,**k)
        for train_ndx, test_ndx in cv.split(X=X_in, y=y_in, groups = groups):
            X_train = X_in[train_ndx, :]
            y_train = y_in[train_ndx]
            w_train = w_in[train_ndx]
            y_test = y_in[test_ndx]

            clf.fit(X=X_train,y=y_train,sample_weight=w_train)
            y_hat = clf.predict_proba(X=X_in[test_ndx, :])
            y_hat = y_hat[::4,]
            w_test = w_in[test_ndx]
            w_test = w_test.reshape((len(w_test)/4,4))
            score = np.sum((y_hat-w_test)**2)
            results = out_results.get(i, [])
            results.append(score)
            out_results.update({i: results})

    for k, v in out_results.items():
        mean_score = np.mean(v)
        out_results.update({k: mean_score})

    best_score = min(out_results.values())
    best_param = min(out_results, key=out_results.get)
    
    return best_score, best_param


def trainRF_CV(pwms, edges, uniqueProteins, obsGrps, full_X, grpInd,
               default = 0):
    # Performs cross-validation and uses either grid searches or random search 
    # to find best set of rfs based on grouped 5-fold cross-validation
    # If default == 1, then just does not perform hyperparamerter tuning

    start, rev = {}, {}
    for k in uniqueProteins:
        start[k], rev[k] = 0, 0

    full_y = formGLM_Y(uniqueProteins)
    full_W = formGLM_trainW(pwms, uniqueProteins, start, rev)

    # For hyperparamter search, create folds where all proteins with
    # identical core sequence are in the same group.
    groups = []
    for i, (s, e) in enumerate(sorted(grpInd.values())):
        for j in range(s,e+1):
            groups.append(i)
    gkf = GroupKFold(n_splits = 5)
    gkf.get_n_splits(full_X[0],full_y[0],groups)
    #rs=RandomizedSearchCV(clf,parameters,scoring='roc_auc',cv=gkf,n_iter=10)    
    
    # Set params to search
    rf_params = {}
    if default:
        for j in range(MWID):
            rf_params[j] = {'max_depth': [None]}
    else:
        max_depth = [3,6,9,12,15,18,None]
        max_features = {}
        for j in range(MWID):
            nf = 20*len(edges[j])
            max_features[j] = [int(round(x*nf)) for x in (0.2,0.3,0.4,0.5,0.6)]
            rf_params[j] = {'max_features': max_features[j],
                            'max_depth': max_depth}
    
    cv_best = {}
    cv_model = {}
    for j in range(MWID):
        param_grid = ParameterGrid(rf_params[j])
        #print param_grid
        fit_params = {'sample_weight': full_W[j]}
        print "Performing grid search for position %d" %j
        if USE_MSE_CV:
            cv_best[j] = \
                grid_cv_mse(full_X[j],full_y[j],full_W[j],gkf,param_grid,
                            groups = groups)
        else:
            cv_best[j] = \
                grid_cv(full_X[j],full_y[j],full_W[j],gkf,param_grid,
                        groups = groups, use_weighting = True)
        #print cv_best[j]
        print list(param_grid)[cv_best[j][1]]
        cv_model[j] = \
            RandomForestClassifier(n_estimators = 100, criterion="entropy",
                                   warm_start=False,n_jobs=-1,
                                   **list(param_grid)[cv_best[j][1]])
        cv_model[j].fit(X=full_X[j],y=full_y[j],sample_weight=full_W[j])
        """
        cv_model[j] = \
            GridSearchCV(RandomForestClassifier(),rf_params[j],
                         cv=gkf,return_train_score = True,
                         scoring = 'neg_log_loss')
        cv_model[j].fit(full_X[j],full_y[j],groups=groups,
                        sample_weight = full_W[j])
        """
    return cv_model

def makePCCtable2(exp, pred, fname, nbrWts = None, 
                  transCols = None, aaPosList = None):
    # Compute the per-position pearson correlation coefficients
    # and mean squared errors between the predicted and the 
    # experimental pwms assuming that the alignment is correct

    fout = open(fname, 'w')
    fout.write('prot\tpos\tpcc\trmse\tic.exp\tic.pred\tnnbrs')
    if transCols is not None:
        fout.write('\ttransfer\tvarPos\n')
    else:
        fout.write('\n')
    for prot in sorted(exp.keys()):
        for pos in range(len(exp[prot])):
            x, y = exp[prot][pos], pred[prot][pos]
            pcc = scipy.stats.pearsonr(x,y)[0]
            rmse = np.sqrt(np.mean((x - y)**2))
            icExp, icPred = information_content(x), information_content(y)
            fout.write('%s\t%d\t%f\t%f\t%f\t%f' %(prot, pos+1, pcc, rmse, icExp, icPred))

            if nbrWts is not None and nbrWts.has_key(prot) and nbrWts[prot].has_key(pos):
                nnbrs = len(nbrWts[prot][pos].keys())
            else:
                nnbrs = 0

            if transCols is None:
                fout.write('\t%d\n' %(nnbrs))
            else:
                if pos in transCols[prot]['wtCols']:
                    transfer = True
                    varPos = aaPosList[transCols[prot]['varPos']]
                else:
                    transfer = False
                    varPos = aaPosList[transCols[prot]['varPos']]
                fout.write('\t%d\t%s\t%d\n' %(nnbrs,transfer,varPos))
    fout.close()

def main():

    np.random.seed(RAND_SEED)

    # Obtain information about the PWM alignment
    outLabel = 'cisbp-chu/structFixed1_grpHoldout_multinom_chain100maxIter15'
    filename = '../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLEFalseChain100Iter15.pickle'
    with open(filename) as f:
        res = pickle.load(f)
    score = [x['ll'] for x in res]
    opt = np.argmax(score)
    reorient = [x['reorient'] for x in res][opt]
    start = [x['start'] for x in res][opt]
    rev = [x['rev'] for x in res][opt]

    # Get the oritented and aligned 6bp pwms from the training set
    # and information about the edges in the model and positions 
    # in the core sequences
    pwms, core, full, edges, edges_hmmPos, aaPosList, testProts = \
        getTrainPairsAndInfo(rescalePWMs = True)
    print len(pwms), len(start)
    assert len(pwms) == len(start)
    pwms_oriented = getOrientedPWMs(pwms, rev)
    pwms_aligned = getAlignedPWMs(pwms_oriented,core,start,MWID,flipAli=False)
    pwms = pwms_aligned

    # Switch to the MIp ranked features from Christensen et al., NAR, 2012
    #for k in core.keys():
    #    core[k] = ''.join([full[k][i] for i in MIP_CORE_POS])
    #print core

    # If want to train using all values in the core sequence for every base position
    if not USE_STRUCT_EDGES:
        allCorePos = set()
        for k in edges.keys():
            allCorePos.add(k)
        for k in edges.keys():
            edges[k] = sorted(list(allCorePos))
    #print edges

    # Set the output directory
    dir = "../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLE"+str(INIT_ORACLE)+"Chain"+\
        str(N_CHAINS)+"Iter"+str(MAXITER)+'/randomForest_scaled50'
    if USE_STRUCT_EDGES:
        dir += '/structEdges'
    else:
        dir += '/allEdges'
    if DEFAULT_PARAMS:
        dir += '_default'
    if USE_MSE_CV:
        dir += '_mseCV'
    if not os.path.exists(dir):
        os.makedirs(dir)

    print("number of proteins used", len(pwms.keys()))
    # Assign cores to similarity groups
    obsGrps = assignObsGrps(core, by = OBS_GRPS)
    with open(dir+'/obsGrpTab.txt', 'w') as fout:
        for k in sorted(obsGrps.keys()):
            for x in obsGrps[k]:
                fout.write('%s\t%s\n' %(k,x))

    print("Output written to: %s" %dir)
    # Write the PWMs used and other basic info to files
    makePWMtab(pwms, dir+'/pwmTab.txt')
    with open(dir+'/coreSeqTab.txt', 'w') as fout:
        fout.write('\n'.join(['%s\t%s' %(k,core[k]) \
                             for k in sorted(core.keys())]))
    
    startTime = time.time()
    uniqueProteins = pwms.keys()
    fullX, grpInd = formGLM_fullX(core,edges,uniqueProteins,obsGrps,numAAs=20)

    #Sanity checks and setting up correct fixed unique protein ordering 
    for g in obsGrps.keys():
        assert len(obsGrps[g])*4 == grpInd[g][1]-grpInd[g][0]+1
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    assert len(uprots) == len(uniqueProteins)
    # So that uniqueProteins is in the same order as obsGrps keys
    uniqueProteins = uprots

    # Train the classifier and pickle the results
    picklefile = dir+'.pickle'
    if FIT_MODEL:
        cv_model = trainRF_CV(pwms, edges, uniqueProteins, obsGrps, fullX, 
                              grpInd, default = DEFAULT_PARAMS)
        print("Ran in %.2f seconds" %(time.time()-startTime))
        #print("Writing results in ", dir+'.pickle')
        #with open(picklefile, 'wb') as f:
        #    pickle.dump(cv_model, f)

    else:
        with open(picklefile) as f:
            clf_res = pickle.load(f)
            cv_model = clf_res

    #print "Best model was selected using CV with scoring as:", cv_model[0].scorer_

    # Extract the best estimator according to the CV criterion
    #print cv_model[0].cv_results_['split0_train_score']
    #print cv_model[0].cv_results_['split0_train_score']
    #print cv_model[0].cv_results_['split1_train_score']
    #print cv_model[0].cv_results_['split1_train_score']
    #print cv_model[0].best_params_


    #"""
    model = {}
    for j in cv_model.keys():
        model[j] = cv_model[j]
        #model[j] = cv_model[j].best_estimator_
        #print "model %d:" %j
        #print cv_model[j].best_params_

    # Make predicitons for the training set and measure performance
    pred_pwms = {}
    for idx, p in enumerate(uniqueProteins):
        pwm = []
        testX = formGLM_testX(fullX, idx)
        for j in range(MWID):
            prediction = model[j].predict_proba(testX[j])[0].tolist()
            pwm.append(prediction)
        pred_pwms[p] = np.array(pwm)
    print("Creating predicted logos according to model ...")
    if MAKE_LOGOS:
        logoDir = dir+'/predicted_logos/'
        if not os.path.exists(logoDir):
            os.makedirs(logoDir)
        makeAllLogos(pred_pwms, core, logoDir)

    # Output the fit information for individual PWMs under S
    pccTabfile = dir+'/pccTable_underS.txt'
    makePCCtable(pwms, pred_pwms, pccTabfile)

    # Make predicted PWMs based on each of the test proteins

    # Get the ground truth test proteins and PWMs
    # Remove examples from cis-bp and chu that correspond to 
    # core seqs in the test set
    fin = open('../hd1-explore/0_splitChu2012-trainTest/testProts.txt','r')
    testProts = [x.strip() for x in fin.readlines()]
    fin.close()
    pwms_test, core_test, full_test = getTestProtsAndPWMs(testProts)
    # Switch to the MIp ranked features from Christensen et al., NAR, 2012
    #for k in core_test.keys():
    #    core_test[k] = ''.join([full_test[k][i] for i in MIP_CORE_POS])
    obsGrps = assignObsGrps(core_test, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots 
    fullX, grpInd = formGLM_fullX(core_test, edges, uniqueProteins, 
                                  obsGrps, numAAs=20)
    pred_pwms = {}
    for idx, p in enumerate(uniqueProteins):
        pwm = []
        testX = formGLM_testX(fullX, idx)
        for j in range(MWID):
            prediction = model[j].predict_proba(testX[j])[0].tolist()
            pwm.append(prediction)
        pred_pwms[p] = np.array(pwm)

    pwms_test_pred = {}  # test pwms aligned to predicted pwms
    for p in pred_pwms.keys():
        # Align to the predicted PWM
        score, shift, ori = comp_matrices(pred_pwms[p], pwms_test[p],
                                          minWidthM2 = 6, oneSided=True)
        shift = -shift
        if ori:
            pwms_test_pred[p] = matrix_compl(pwms_test[p])[shift:shift+MWID]
        else:
            pwms_test_pred[p] = pwms_test[p][shift:shift+MWID]
    outDir = dir+'/testFit/'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    makePCCtable2(pwms_test_pred, pred_pwms, outDir+'pccTable_test_pred.txt')
    #"""

if __name__ == '__main__':
    main()


