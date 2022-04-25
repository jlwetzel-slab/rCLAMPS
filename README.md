# rCLAMPS

rCLAMPS is software implementing the method described in the manuscript:  "Learning probabilistic protein-DNA recognition codes from DNA-binding specificities using structural mappings" by Joshua Wetzel, Kaiqian Zhang, and Mona Singh.

The rCLAMPS framework implements a Gibbs sampling approach to examine protein-DNA interaction data (in the form of paired position-specific weight matrices (PWMs) and corresponding protein sequences) for a DNA-binding family of proteins, along with aggregated structural information for the that protein family.  It simultaneously infers an interpretable predictive model of the protein family's amino acid-to-base interaction preferences as well as a set of inferred corresponding protein-DNA interfaces via alignment to a canonical amino acid-to-base position "contact map" for the protein faimily.

A short script reproducing analysis plots and and other data reported in the manuscript is located in ./code/analysis_manuscript.R.  The in-code documentation there describes which additional code files were run to produce inputs to that script.

The main file of interest for running the rCLAMPS framework itself is ./code/gibbsAlign_GLM.py - it provides the following input/output functionalities:

Inputs:
1.  A protein-DNA interface contact map for the protein family of interest - I.e., positions within the DNA-binding domain that canoncially interact with DNA bases, along with the corresponding DNA position indices that they tend to interact with (according to a canonical numbering scheme, based on the underlying multiple structural alignment).
2.  A set of protein-DNA binding specificities as PWMs, along with their corresponding protein sequences as a fasta file.

To reproduce the model described in the manuscript above, the software computes these inputs using various data files stored in this repo, but the code also accepts tab-delimited text file inputs (see inputs in the ./precomputedInputs/ directory) by setting global variables appropriately at the top of the script.  Please see the ./code/gibbsAlign_GLM.py code documentation for details.
 
./code/gibbsAlign_GLM.py outputs a pickle file containing a dictionary of Python list objects (each index in the list corresponds to one of K Gibbs sampling chains), keyed as follows, along with a set of predicted motifs based on the optimal 'final_model' object (see below):
1.  'final_model':  A dictionary of scikit-learn LogisticRegression objects with multiclass='multinomial', one for each base position in the contact map, keyed by the base position.  Auxilliary functions needed to transform protein sequence inputs into the proper input format for the models and to make specificity predicitons for novel proteins are included in gibbsAlign_GLM.py.
2.  'll':  The log likelihood of the final model.
3.  'start':  A dictionary of start positions, keyed by the PWM name, of the protein-DNA interaction interface inferred by while estimating the final_model, assuing the PWMs are oriented in the direction given by 'rev'.
4.  'rev':  A dictionary of boolean values, keyed by the PWM name, of the PWM orientations inferred by while estimating the final_model.  0 is the original orientation, 1 is the reverse complement orientation.
5.

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

Prerequisite software and packages:
1.  Python 2.7 with the scikit-learn and scipy packages
2.  HMMer3, for mapping proteins sequences to domain positions in the protein-DNA structural interaction interface model.  The code has been tested up to version HMMer3.3.1, though it should work with any version provided the output formatting has not changed.

Currently the software is tested on a Linux machine running Python v.2.7.17 with scikit-learn v.0.20.4, numpy v1.16.6, and scipy v.1.2.3.  In order to run the code to reproduce the model from the manuscript, install the above software on your local machine, git clone this directory, and set the following global variable at the top of ./code/gibbsAlign_GLM.py: HMMER_HOME = 'path to the hmmsearch executable from HMMer3'. Then simply run from the ./code subdirectory:  python gibbsAlign_GLM.py

This will output to the 'new_output' subdirectory.  Note that as currently set up, this will run 100 gibbs sampling chains in parallel (using up to 1 - number of CPUs on your machine).  This may take some time depending on the number of CPUs/cores on your machine.  The number of sampling chains and maximum number of iterations per chain can be changed using the global variables N_CHAINS and MAXITER, respectively.
