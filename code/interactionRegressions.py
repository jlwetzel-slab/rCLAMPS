# Make regression models with interaction terms based on the PWMs that were 
# alligned by the training procedure

from gibbsAlign_GLM import *
from hd1_transfer_predictions import getTrainPairsAndInfo, getTestProtsAndPWMs
from sklearn.ensemble import RandomForestRegressor
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
MAKE_LOGOS = False

################

HMM_OFFSET = 2
MIP_CORE_POS = [3,6,19,47,50,54,55]
MIP_CORE_POS = [x - HMM_OFFSET for x in MIP_CORE_POS]
FIT_MODEL = False
USE_STRUCT_EDGES = False
DEFAULT_PARAMS = False
USE_MSE_CV = False
RESCALE_PWMS = True
PERFORM_NESTED_CV = False

def trainRF_CV(pwms, edges, uniqueProteins, obsGrps, full_X, grpInd,
               default = 0):
    # Performs cross-validation and uses either grid searches or random search 
    # to find best set of rfs based on grouped 5-fold cross-validation
    # If default == 1, then just does not perform hyperparamerter tuning

    start, rev = {}, {}
    for k in uniqueProteins:
        start[k], rev[k] = 0, 0   
    full_y = formGLM_trainW(pwms, uniqueProteins, start, rev, 
                            modelType='regressor')
    #print(full_y[2]['A'])
    #print(full_y[0]['A'])

    # For hyperparamter search, create folds where all proteins with
    # identical core sequence are in the same group.
    groups = []
    for i, (s, e) in enumerate(sorted(grpInd.values())):
        for j in range(s,e+1):
            groups.append(i)
    gkf = GroupKFold(n_splits = 5)
    gkf.get_n_splits(full_X[0],full_y[0]['A'],groups)
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
            max_features[j] = [int(round(x*nf)) for x in (0.2,0.3,0.4,0.5,0.6)]+["auto"]
            rf_params[j] = {'max_features': max_features[j],
                            'max_depth': max_depth}  
    cv_best = {}
    cv_model = {}
    for j in range(MWID):
        cv_model[j] = {}
        for base in BASE:
            print "Performing grid search for position %d, base %s" %(j,base)
            cv_model[j][base] = \
                GridSearchCV(RandomForestRegressor(n_estimators=50),
                             param_grid=rf_params[j],cv=gkf,return_train_score=True,
                             refit=True,iid=True)
            cv_model[j][base].fit(full_X[j],full_y[j][base],groups=groups)
    return cv_model

def trainRF_nestedCV(pwms, core, edges, uniqueProteins, obsGrps, full_X, grpInd, outDir):
    # Performs cross-validation and uses either grid searches or random search 
    # to find best set of rfs based on grouped 5-fold cross-validation
    # If default == 1, then just does not perform hyperparamerter tuning

    start, rev = {}, {}
    for k in uniqueProteins:
        start[k], rev[k] = 0, 0   
    full_y = formGLM_trainW(pwms, uniqueProteins, start, rev, 
                            modelType='regressor')

    # For hyperparamter search, create folds where all proteins with
    # identical core sequence are in the same group.
    groups = []
    for i, (s, e) in enumerate(sorted(grpInd.values())):
        for j in range(s,e+1):
            groups.append(i)  
    groups = np.array(groups)

    #"""   
    # Set up the params to search
    rf_params = {}
    max_depth = [3,6,9,12,15,18,None]
    max_features = {}
    for j in range(MWID):
        nf = 20*len(edges[j])
        max_features[j] = [int(round(x*nf)) for x in (0.2,0.3,0.4,0.5,0.6)]+["auto"]
        rf_params[j] = {'max_features': max_features[j],
                        'max_depth': max_depth}  

    # Get the nested k-fold CV performance
    gkf_outer = GroupKFold(n_splits = 10)
    gkf_outer.get_n_splits(full_X[0],full_y[0]['A'],groups)
    cv_model = {}
    k = 0
    for train_ndx, test_ndx in gkf_outer.split(full_X[0], groups = groups):
        print "Training optimal model for outer CV fold %d ..." %k

        # Subset the protein list, core dicts, and pwm dicts for this split
        test_cores, test_pwms = {}, {}
        for i, prot in enumerate(uniqueProteins):
            if i in test_ndx:
                test_cores[prot] = core[prot]
                test_pwms[prot] = pwms[prot]

        # Train the model for this split
        model = {}
        for j in range(MWID):
            # split data
            model[j] = {}
            X_train, X_test = full_X[j][train_ndx,:], full_X[j][test_ndx,:]
            grps_inner = groups[train_ndx]
            gkf_inner = GroupKFold(n_splits = 5)
            gkf_inner.get_n_splits(X_train, groups = grps_inner)
            for base in BASE:
                y_train,y_test = full_y[j][base][train_ndx],full_y[j][base][test_ndx]
                print "Performing grid search for position %d, base %s" %(j,base)
                search = \
                    GridSearchCV(RandomForestRegressor(n_estimators=50),
                                 param_grid=rf_params[j],cv=gkf_inner,
                                 refit=True,iid=True)
                res = search.fit(X_train,y_train,groups=grps_inner)
                model[j][base] = res.best_estimator_
        
        # Make predicitons for the held out test fold
        pred_pwms, pwms_test_pred = \
            makePredsAndAlign(test_cores, test_pwms, model, edges)

        # Analyze performance on the held-out fold in terms of PCC and RMSE
        makePCCtable2(pwms_test_pred, pred_pwms, 
                      outDir+'pccTable_test_pred_fold%d.txt' %k)
        k += 1

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

def makeRegressorPreds(models, X, protLabs):
    # Returns the combined predictions for the regression models
    # given amino acid input vectors in dictionary X, keyed per 
    # model position

    pred_pwms = {}
    for bpos in range(MWID):
        pred_pwms[bpos] = np.zeros((len(X[bpos]),4), dtype = float)
        for base in BASE:
            preds = models[bpos][base].predict(X[bpos])
            pred_pwms[bpos][:,B2IND[base]] = preds
        pred_pwms[bpos] = pred_pwms[bpos]/pred_pwms[bpos].sum(axis=1,keepdims=True)
        #print bpos, pred_pwms[bpos]

    perProt = {}
    for i, p in enumerate(protLabs):
        perProt[p] = np.zeros((MWID,4), dtype = float)
        for bpos in pred_pwms.keys():
            perProt[p][bpos,:] = pred_pwms[bpos][i]
    #print perProt
    return perProt

def make_testX(fullX, uniqueProteins):
    # Returns a dictionary protein encodings
    protVects = {}
    testX = {}
    for idx, p in enumerate(uniqueProteins):
        protVects[p] = formGLM_testX(fullX, idx, modelType = 'regressor')
    for bpos in range(MWID):
        testX[bpos] = np.stack([protVects[p][bpos][0] for p in uniqueProteins])
    return testX

def makePredsAndAlign(core_test, pwms_test, model, edges):
    # Makes PWM predictions a test set based on the set of RF models (model), 
    # aligns to the ground truth PWMs, and returns the predicted PWMs
    # and the aligned (and truncated) ground truth PWMs

    # Switch to the MIp ranked features from Christensen et al., NAR, 2012
    #for k in core_test.keys():
    #    core_test[k] = ''.join([full_test[k][i] for i in MIP_CORE_POS])
    obsGrps = assignObsGrps(core_test, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots 
    fullX, grpInd = formGLM_fullX(core_test, edges, uniqueProteins, 
                                  obsGrps, numAAs=20, modelType='regressor')

    # Make predicitons for the test set and measure performance
    testX = make_testX(fullX, uniqueProteins)
    pred_pwms = makeRegressorPreds(model, testX, uniqueProteins)
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
    return pred_pwms, pwms_test_pred

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
        getTrainPairsAndInfo(rescalePWMs = RESCALE_PWMS)
    print len(pwms), len(start)
    assert len(pwms) == len(start)
    pwms_oriented = getOrientedPWMs(pwms, rev)
    pwms_aligned = getAlignedPWMs(pwms_oriented,core,start,MWID,flipAli=False)
    pwms = pwms_aligned

    #"""
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

    # Set the output directory
    dir = "../results/cisbp-chu/structFixed1_grpHoldout_multinomial_ORACLE"+str(INIT_ORACLE)+"Chain"+\
        str(N_CHAINS)+"Iter"+str(MAXITER)+'scaled50/randomForest_regression_scaled50'
    if USE_STRUCT_EDGES:
        dir += '/structEdges'
    else:
        dir += '/allEdges'
    if DEFAULT_PARAMS:
        dir += '_default'
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
    fullX, grpInd = formGLM_fullX(core,edges,uniqueProteins,obsGrps,
                                  numAAs=20,modelType='regressor')
    print fullX[0].shape, len(grpInd.keys())

    #Sanity checks and setting up correct fixed unique protein ordering 
    for g in obsGrps.keys():
        assert len(obsGrps[g]) == grpInd[g][1]-grpInd[g][0]+1
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    assert len(uprots) == len(uniqueProteins)
    # So that uniqueProteins is in the same order as obsGrps keys
    uniqueProteins = uprots

    # Do nested CV so that we can compare performance across
    # model types independently of performance on the test set
    if PERFORM_NESTED_CV:
        outDir = dir+'/nestedCV_fit/'
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        trainRF_nestedCV(pwms, core, edges, uniqueProteins, obsGrps, fullX, 
                         grpInd, outDir)

    # Train the classifier, finding the hyperparameters with best average
    # performance across 5 test folds, and save the resulting model
    picklefile = dir+'.pickle'
    if FIT_MODEL:
        cv_model = trainRF_CV(pwms, edges, uniqueProteins, obsGrps, fullX, 
                              grpInd, default = DEFAULT_PARAMS)
        print("Ran in %.2f seconds" %(time.time()-startTime))
        print("Writing results in ", dir+'.pickle')
        with open(picklefile, 'wb') as f:
            pickle.dump(cv_model, f)

    else:
        with open(picklefile) as f:
            clf_res = pickle.load(f)
            cv_model = clf_res

    # Extract the models with highest average CV performance across folds
    model = {}
    for j in cv_model.keys():
        model[j] = {}
        for base in BASE:
            model[j][base] = cv_model[j][base].best_estimator_
            print cv_model[j][base].best_params_

    # Make predicitons for the training set and measure performance
    testX = make_testX(fullX, uniqueProteins)
    pred_pwms = makeRegressorPreds(model, testX, uniqueProteins)
    #print pred_pwms

    # Output the fit information for individual PWMs under S
    pccTabfile = dir+'/pccTable_underS.txt'
    makePCCtable(pwms, pred_pwms, core, pccTabfile)

    # Make predicted PWMs based on each of the test proteins
    with open('../hd1-explore/0_splitChu2012-trainTest/testProts.txt','r') as fin:
        testProts = [x.strip() for x in fin.readlines()]
    pwms_test, core_test, full_test = getTestProtsAndPWMs(testProts)
    pred_pwms, pwms_test_pred = makePredsAndAlign(core_test, pwms_test, model, edges)


    outDir = dir+'/testFit/'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    makePCCtable2(pwms_test_pred, pred_pwms, outDir+'pccTable_test_pred.txt')
    makePWMtab(pred_pwms,outDir+'/pwms_pred.txt')
    makePWMtab(pwms_test_pred,outDir+'/pwms_testAli.txt')

    if MAKE_LOGOS:
        logoDir = dir+'/predicted_logos/'
        print("Creating predicted logos according to model ...")
        if not os.path.exists(logoDir):
            os.makedirs(logoDir)
        makeAllLogos(pred_pwms, core, logoDir)

    #"""

if __name__ == '__main__':
    main()


