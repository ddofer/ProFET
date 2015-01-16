To use the entire pipeline from the command line, enter:
"python pipeline.py" from the directory containing "pipeline.py" with the appropriate parameters.

The parameters are:
--trainingSetDir					The path to the training set fasta files
usage example: python pipeline.py --trainingSetDir r'C:\Feature_Extract\test_seq\Chap\train'
--testingDir						The path to the testing set fasta files
usage example: python pipeline.py --testingDir r'C:\Feature_Extract\test_seq\Chap\test'
--resultsDir						The path to the directory where the results will be saved
usage example: python pipeline.py --resultsDir r'C:\Feature_Extract\test_seq\Chap\results'
--trainFeatures						A flag indicating whether to extract (and save) the training set features
usage example: python pipeline.py --trainFeatures True
--testFeatures						A flag indicating whether to extract (and save) the testing set features
usage example: python pipeline.py --testFeatures False
--classType							The type of class the protein will be sorted by. The three options are:
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
python pipeline.py --trainingSetDir r'C:\Feature_Extract\test_seq\Chap\train' --testingDir r'C:\Feature_Extract\test_seq\Chap\test' --trainFeatures True --testFeatures True --classType 'dir'