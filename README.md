# rCLAMPS

rCLAMPS is software implementing the method described in the manuscript:  "Learning probabilistic protein-DNA recognition codes from DNA-binding specificities using structural mappings" by Joshua Wetzel, Kaiqian Zhang, and Mona Singh.

The rCLAMPS framework implements a Gibbs sampling approach to examine protein-DNA interaction data (in the form of paired position-specific weight matrices (PWMs) and corresponding protein sequences) for a DNA-binding family of proteins, along with aggregated structural information for that protein family.  It simultaneously infers an interpretable predictive model of the protein family's amino acid-to-base interaction preferences (recognition code) as well as a set of corresponding protein-DNA interfaces relative to a canonical amino acid-to-base position "contact map" for the protein faimily.  We have tested the software on the homoedomain proteins and a subset of the C2H2-ZFs for which it is known in advance which domains with the C2H2-ZF protein contact DNA.

If you are interested in making predictions for novel homoedomain proteins, we provide an example for using the homeomdomain recognition code inferred by rCLAMPS in the file ./code/examplePredictions.py.  All that is required for doing so is that you set the global variable FASTA_INPUT_FILE to point to a fasta file wherein each entry is a protein sequence containing a single homeodomain instance (see, e.g.,  ./examplePredictions/predictionExamples.fa) and you set OUTPUT_DIR to point to a desired output directory path.  Note that examplePredictions.py requires that each protein sequence in this fasta file have a unique name, and that each entry contains a single homeodomain instance with all base-contacting HMM match states present.  For each protein in the fasta file that meets these criteria, a predicted PWM in tabular format will be written into a file called OUTPUT_DIR/predicted_pwms.txt.  If you have weblogo installed (see http://weblogo.threeplusone.com/manual.html#download), then you can set MAKE_LOGOS = True and each predicted logo will also be visualized and saved into OUTPUT_DIR/predicted_logos/.

For a list of requirements to run rCLAMPS and/or examplePredictions.py, please see requirements.txt. If yu wish to reproduce models examined in the manuscript, please see modelReproductions.txt.  A short script reproducing analysis plots reported in our main manuscript are located in ./code/analysis_manuscript.R. and abalysis_manuscript_zf-C2H2.R The in-code documentation there describes which additional code files were run to produce inputs the tables analysed in that script.

If you wish to run the rCLAMPS framework itself to infer new recogntion codes or train with different datasets, the main file of interest is ./code/gibbsAlign_GLM.py - it provides the following input/output functionalities:

Input types:
1.  A protein-DNA interface contact map for the protein family of interest - I.e., positions within the DNA-binding domain that canoncially interact with DNA bases, along with the corresponding DNA position indices that they tend to interact with (according to a canonical numbering scheme, based on the underlying multiple structural alignment).  See ./precomputedInputs/homeodomain_contactMap.txt for an example.
2.  A set of protein-DNA binding specificities as PWMs, along with their corresponding protein sequences as a fasta file.  See ./precomputedInputs/pwmTab_homeodomains_all.txt and ./precomputedInputs/proteins_homeodomains_hasPWM.fa for examples. 
 
./code/gibbsAlign_GLM.py outputs a pickle file containing a dictionary of Python list objects (each index in the list corresponds to one of K Gibbs sampling chains), keyed as follows, along with a set of predicted motifs based on the optimal 'final_model' object (see below):
1.  'final_model':  A dictionary of scikit-learn LogisticRegression objects with multiclass='multinomial', one for each base position in the contact map, keyed by the base position.  Auxilliary functions needed to transform protein sequence inputs into the proper input format for the models and to make specificity predicitons for novel proteins are included in gibbsAlign_GLM.py.
2.  'll':  The log likelihood of the final model.
3.  'start':  A dictionary of start positions, keyed by the PWM name, of the protein-DNA interaction interface inferred by while estimating the final_model, assuing the PWMs are oriented in the direction given by 'rev'.
4.  'rev':  A dictionary of boolean values, keyed by the PWM name, of the PWM orientations inferred by while estimating the final_model.  0 is the original orientation, 1 is the reverse complement orientation.

For example, to extract the set of optimal models, starts, and orientations from the pickle file, one would use: 

```
with open(filename) as f:
 res = pickle.load(f)
score = [x['ll'] for x in res]
opt = np.argmax(score)
start = [x['start'] for x in res][opt]
rev = [x['rev'] for x in res][opt]
opt_models = [x['final_model'] for x in res][opt]
```

This code is feely available for reuse and modification.  Please refer any questions about the software to jlwetzel@princeton.edu and/or mona@cs.princeton.edu.
