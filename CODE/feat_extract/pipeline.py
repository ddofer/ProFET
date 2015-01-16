#!/sw/bin/python3.3

__author__ = 'Michael+Dan'
'Last Modified from Michael doron - 25.11.2014'
'''
TODO: Use OS.path,join insteado f + strings. (cross OS compatability + saves headaches.
Also - make it clearer about how to just predict for example.. and valid/invalid options
'TODO: Add "Get top features"  to command line supported options. And GetTraining perf (via pipeline tasks)'

TODO: Add option for user to choose what model to use for predictions! (And params..)
'''


from FeatureGen import featExt
from Model_trainer import trainClassifier
from FeatureGen import writeClassifiedFastas

import time
import pandas as pd
import sklearn
import numpy
from sklearn.preprocessing import MinMaxScaler, StandardScaler, LabelEncoder

import click


@click.command()
@click.option('--trainingSetDir','-r','trainingDir',default=r'C:\Users\Michael\Desktop\Feature_Extract\test_seq\Chap\train',
              help='The path to the training set fasta files', type = str)
@click.option('--testingSetDir','-s','testingDir',default=r'C:\Users\Michael\Desktop\Feature_Extract\test_seq\Chap\test',
              help='The path to the testing set fasta files', type = str)
@click.option('--resultsDir','-rs','resultsDir',default=r'C:\Users\Michael\Desktop\Feature_Extract\test_seq\Chap\results',
              help='The path to directory to write the results files', type = str)
@click.option('--trainFeatures','-rf','GetTrainingFeatures',default=True,help='Whether to extract the training set features', type = bool)
@click.option('--testFeatures','-sf','GetTestFeatures',default=False,help='Whether to get the testing set features', type = bool)
@click.option('--classType','-ct','classType',default='file',help="Defines the classname of each protein, by \'dir\', \'file\', or \'id\'.", type = str)
def pipeline(trainingDir,testingDir,resultsDir, GetTrainingFeatures,GetTestFeatures, classType):
    # change here to the training data folder
    # trainingDir = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\CODE\Feature_Extract\test_seq\Chap'
    # change here to the testing data folder
    # testingDir = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\FASTA_Sets\HSP33_Chap\Unknown_Tests'

    if GetTrainingFeatures==True:
        print('Starting to extract features from training set')
        'Temporary measure: If features extracted and saved, disable following line to avoid re-extracting trainign features'
        featExt(directory=trainingDir, trainingSetFlag=True,
                classType=classType, normParams='.')
        print('Extracted training data features')


    'TODO: Seperate model training/prediction from feat.extraction!'
    if GetTestFeatures==True:
        print('Training predictive model')
        #ORIG \\
##        model, lb_encoder = trainClassifier(trainingDir+'\\trainingSetFeatures.csv', False, 'forest', 0, False, False)
        model, lb_encoder = trainClassifier(filename=trainingDir+'/trainingSetFeatures.csv',normFlag= False,classifierType= 'forest',kbest= 0,alpha= False,optimalFlag= False)
        print('Model trained')

    'Change to "If GetPredictions==True" , after adding such a param'
    if GetTestFeatures==True:
            print('Extracting features from test set')
            featExt(testingDir, False, 'dir', trainingDir+'\\trainingSetNormParams.csv')
            print('Extracted test data features')
            dfTesting = pd.DataFrame.from_csv(testingDir+'\\testingSetFeatures.csv')
            dfTraining = pd.io.parsers.read_csv(trainingDir+'/trainingSetFeatures.csv',nrows=2)
            '''
            # FeatureFilt
            Filter Extracted Features, keeping only feats that are in the training set.
            This is crucial! (Remember  be reapplied if used elsewhere, if feature filtering/selection used)
            '''
            # remove feature in dfTesting when not in dfTesting
            feature_cols = [col for col in dfTraining.columns if col not in ['classname','Id','proteinname']]
            features = dfTesting[feature_cols].values

            print('Predicting labels')
            results = model.predict(features)
            labels = lb_encoder.inverse_transform(results)
            # dfTesting['classname'].append(list(labels))
            dfTesting['classname'] = labels
            #df to df2 :
            df2 = dfTesting['classname']
            df2.to_csv(testingDir+'\\PredictedTestSetResults.csv')
            print('saved results to ' + testingDir+'\\PredictedTestSetResults.csv')
            writeClassifiedFastas(classType, testingDir, resultsDir, df2)
  
if __name__ == '__main__' :
    pipeline()
