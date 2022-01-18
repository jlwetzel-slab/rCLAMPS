# rCLAMPS

rCLAMPS is software implementing the method described in the manuscript:  "Learning probabilistic protein-DNA recognition codes from DNA-binding specificities using structural mappings" by Joshua Wetzel, Kaiqian Zhang, and Mona Singh.

The software implements a Gibbs sampling approach to examine protein-DNA interaction data (in the form of paired position-specific weight matrices (PWMs) and corresponding protein sequences) for a DNA-binding family of proteins, along with aggregated structural information for the that protein family, to simultaneously infer an interpretable predictive model of the protein family's amino acid-to-base interaction preferences as well as a set of inferred corresponding protein-DNA interfaces.

The main file of interest for running the rCLAMPS method is ./code/gibbsAlignGLM.py - it provides the following input/output functionalities:

Inputs:
1.  A protein-DNA structural interaction interface model for the protein family of interest - I.e., positions within the DNA-binding domain that interact with DNA bases, along with the corresponding DNA position indices that they interact with (according to a canonical numbering scheme).
2.  A set of protein-DNA binding specificities as PWMs, along with their corresponding protein sequences as a fasta file.

To reproduce the model discussed in the manuscript above, the software computes these inputs using various data files stored in this repo, but the code will ultimately be altered to take these inputs modularly from any source in a prescribed format.
 
It outputs a pickle file containing a dictionary of Python list objects (each index in the list corresponds to one of K gibbs sampling chains), keyed as follows:
1.  'final_model':  A dictionary of scikit-learn LogisticRegression objects with multiclass='multinomial', one for each base position in the protein-DNA structural interaction interface model, keyed by the base position.  Auxilliary functions needed to transform protein sequence inputs into the proper input format for the models are included in gibbsAlignGLM.py.
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

Prerequisite software and packages:
1.  Python 2.7 with the scikit-learn and scipy packages
2.  HMMer3, for mapping proteins sequences to domain positions in the protein-DNA structural interaction interface model.  The code has been tested up to version HMMer3.3.1, though it should work with any version provided the output formatting has not changed.

In order to run the code to reproduce the model from the manuscript, install the above software on your local machine, git clone this directory, and set the following global variable at the top of ./code/gibbsAlign_GLM.py: HMMER_HOME = 'path to the hmmsearch executable from HMMer3'. Then simply run from the ./code subdirectory:  python gibbsAlign_GLM.py

This will output to the 'new_output' subdirectory.  Note that as currently set up, this will run 100 gibbs sampling chains in parallel (using up to 1 - number of CPUs on your machine).  The number of sampling chains and maximum number of iterations per chain can be changed using the global variables N_CHAINS and MAXITER, respectively.
