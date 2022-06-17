# A sample file for how to make a prediction for a novel
# homoedomain protein using the rCLAMPS model

from predictionExamples_helpers import *

MODEL_DIR = '../my_results/allHomeodomainProts/'   # Directory containing the rCLAMPS output

# Directory containing a fasta of homeomdomain proteins to predict specificities for
FASTA_INPUT_FILE = '../examplePredictions/predictionExamples.fa'  
OUTPUT_DIR = '../examplePredictions/'

MAKE_LOGOS = True  # Set to True to make logos using WebLogo (False for only text PMWs)

def createFinalModel():

    # Import the training dataset and contact map information in the 
    # format used in '../precomputedInputs/'
    pwms, core, full, edges, edges_hmmPos, aaPosList, testProts = getPrecomputedInputs()
    
    # Assign to distinct observation groups
    obsGrps = assignObsGrps(core, by = OBS_GRPS)
    uprots = []
    for grp in obsGrps.keys():
        uprots += obsGrps[grp]
    uniqueProteins = uprots  

    nDoms = {}
    for p in uniqueProteins:
        nDoms[p] = len(core[p])/len(aaPosList)

    # Train the model using the optimal offset found previously by rCLAMPS
    filename = MODEL_DIR+'result.pickle'
    with open(filename) as f:
        res = pickle.load(f)
    score = [x['ll'] for x in res]
    opt = np.argmax(score)
    start = [x['start'] for x in res][opt]
    rev = [x['rev'] for x in res][opt]

    fullX, grpInd = formGLM_fullX(core, edges, uniqueProteins, obsGrps)
    model = form_model(fullX, uniqueProteins, nDoms, pwms, start, rev)
    return model, aaPosList, edges

def main():

    # Create a temporary directory
    if not os.path.exists('./tmp/'):
        os.makedirs('./tmp')

    # Create the model object and get the list of relevant match states
    print "Creating the model object ..."
    model, aaPosList, edges = createFinalModel()

    print "Converting fasta proteins to X vectors ..."
    # Read in proteins you are interested in predicting specificities for 
    # and create appropriate X matrices
    fullX, uniqueProteins, obsGrps, grpInd = \
        getFullX_fromFasta(FASTA_INPUT_FILE, aaPosList, edges)
    print fullX.keys()

    print "Making PWM predictions ..."
    # Make a prediction for each protein with a distinct combination 
    # of base-contatcing residues and store
    pred_pwms = {}
    for k, coreSeq in enumerate(grpInd.keys()):
        startInd_ho, endIndex_ho = grpInd[coreSeq]
        testProteins = []
        for prot in uniqueProteins:
            if prot in obsGrps[coreSeq]:
                testProteins += [prot]
        testX = formGLM_testX(fullX, startInd_ho, startInd_ho + 4)
        pwm = []
        for j in range(MWID):
            prediction = model[j].predict_proba(testX[j])[0].tolist()
            pwm.append(prediction)
        for p in testProteins:
            pred_pwms[p] = np.array(pwm)

    # Output the specificity predictions in a table format
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    makePWMtab(pred_pwms,OUTPUT_DIR+'/predicted_pwms.txt')

    if MAKE_LOGOS:
        print "Generating logos for each prediciton"
        makeLogos(pred_pwms, OUTPUT_DIR + '/predicted_logos/')

if __name__ == '__main__':
    main()