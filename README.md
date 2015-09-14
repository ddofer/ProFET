# ProFET
ProFET: Protein Feature Engineering Toolkit for Machine Learning


To run ProFET you will need Python 3 and the following packages installed:
* Pandas
* Numpy
* Scikit-learn
* Biopython

I STRONGLY recommend installing them using the Anaconda Python distribution if you don't have them already, It's easiest:

All of these packages are part of the default Anaconda destribution, except for biopython. You can add Biopython by running from the command line:
* conda install Biopython  *

Code for running the package is in the subdirectory 'CODE/feat_extract' .
You can run the feature extraction pipeline using 'pipeline.py'. (Detailed usage instructions can be found in its accompanying readme).

You can also extract features yourself using the underlying methods.
Currently, you can alter the features extracted by 'pipeline.py' by directly modifying the parameters of the function "Get_Protein_Feat" in '/feat_extract/FeatureGen.py'.

The datasets used in the paper are found in the 'fasta' directory. Be careful with file/fasta formats and names.

Please note that the code is currently 'beta' - that means that it's rough, and there will be bugs. Feel free to add or improve!


If you use our tool, or code, please cite us!:
**Ofer, Dan, and Michal Linial. "ProFET: Feature engineering captures high-level protein functions." Bioinformatics (2015): btv345.**

Dan.
