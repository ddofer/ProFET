To use the entire pipeline from the command line, enter:
"python pipeline.py" from the directory containing "pipeline.py" with the appropriate parameters.

The parameters are:
--trainingSetDir					The path to the training set fasta files
usage example: python pipeline.py --trainingSetDir r'C:\Feature_Extract\test_seq\Chap\train'
--testingSetDir						The path to the testing set fasta files
usage example: python pipeline.py --testingDir r'C:\Feature_Extract\test_seq\Chap\test'
--resultsDir						The path to the directory where the results will be saved
usage example: python pipeline.py --resultsDir r'C:\Feature_Extract\test_seq\Chap\results'
--trainFeatures						A flag indicating whether to extract (and save) the training set features
usage example: python pipeline.py --trainFeatures True
--testFeatures						A flag indicating whether to extract (and save) the testing set features
usage example: python pipeline.py --testFeatures True
--classType							The type of class the protein will be sorted by. The 3 options are:
									'dir' - each protein will get the class of its immediate directory name,
										for instanced, a protein whose fasta file is sitting in the directory
										'C:\Feature_Extract\test_seq\Chap\train\A' will get the class 'A'
									'file' - each protein will get the class of its fasta file name.
										for instance, a protein written in the fasta file
										'C:\Feature_Extract\test_seq\Chap\train\A\catProteins.fasta'
										will get the class 'catProteins'
									'id' - each protein will get the class of its header id, for instance
										a protein with the header '>AHG85394.1 Chaperone protein DnaJ [Bibersteinia trehalosi USDA-ARS-USMARC-190]'
										will get the class "Chaperone"
usage example: python pipeline.py --classType 'dir'

To conclude, if you want to train the classifier on fasta files from the directory C:\Feature_Extract\test_seq\Chap\train,
and extract the features from these files (this is necessary for training the classifier),
and test the classifier on fasta files from the directory C:\Feature_Extract\test_seq\Chap\test,
while extracting features from these files (this is necessary for testing the classifier),
and classifying the proteins by their immediate directory name,
you need to enter:

python pipeline.py --trainingSetDir r'C:\Feature_Extract\test_seq\Chap\train' --testingSetDir r'C:\Feature_Extract\test_seq\Chap\test' --trainFeatures True --testFeatures True --classType file

NOTE: The ClassType should be written WITHOUT any ''. (i.e: file , not 'file').
Note - Directory names may not contain any spaces!

Note - Fasta files will be acquired recursively from any subdirectories in the folder!

---------------------------------------------------------------------------------------
Examples:

"Toy" + no unknowns to predict:
python pipeline.py --trainingSetDir /cs/prt3/danofer/ProtFeat/t --classType ‘file’ --testingSetDir /cs/prt3/danofer/ProtFeat/t --trainFeatures True


	To Get a predicted "good" machine learning model and hyperparameters, selected by cross validation, you can use pipetasks.py:
EG your training data has been extracted (and filename is unchanged) to 'cs/prt3/danofer/ProtFeat/t', then you need to enter (at the command line):
python PipeTasks.py /cs/prt3/danofer/ProtFeat/t

For Example:
python pipeline.py --trainingSetDir /cs/prt3/danofer/NP-Collab/npid_2/datasets/final_collab_2014_npp_trainingsets/Train --classType file --trainFeatures True  --resultsDir /cs/prt3/danofer/NP-Collab/npid_2/datasets/final_collab_2014_npp_trainingsets/Train

python PipeTasks.py /cs/prt3/danofer/NP-Collab/npid_2/datasets/final_collab_2014_npp_trainingsets/Train



	To Get the K top features, filtered by alpha, you may use the GetKFeatures function from pipetasks:
EG: From the directory containg the scripts (including PipeTasks.py), at the command line, to use the function on the (previously generated!) file/csv holding the features we wish to reduce (located in "/a/fr-05/vol/protein/danofer/ProtFeat/t" in this case) ; (from in python):

	python
	import PipeTasks
	import os
	os.chdir("/a/fr-05/vol/protein/danofer/ProtFeat/t")
	PipeTasks.GetKFeatures('trainingSetFeatures.csv')