#!/sw/bin/python3.3

__author__ = 'Michael'
'Modified from Michael doron - 4.11.2014'
'''
TODO: Use OS.path,join insteado f + strings. (cross OS compatability + saves headaches.
Also - make it clearer about how to just predict for example.. and valid/invalid options
'TODO: Add "Get top features"  to command line supported options. And GetTraining perf (via pipeline tasks)'

TODO: Add option for user to choose what model to use for predictions!
TODO: Have "getTestSet" also filter features (keep only feats in the training data!)
'''

# from Feature_Extract.FeatureGen import featExt
# from Machine_Learn.Model_trainer import trainClassifier
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
    # trainingDir = str(trainingDir)

    # change here to the testing data folder
    # testingDir = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\FASTA_Sets\HSP33_Chap\Unknown_Tests'
    # testingDir = str(testingDir)

    if GetTrainingFeatures==True:
        print('Starting to extract features from training set')
        'Temporary measure: If features extracted and saved, disable following line to avoid re-extracting trainign features'
        featExt(directory=trainingDir, trainingSetFlag=True,
                classType=classType, normParams='.')
        print('Extracted training data features')


    'TODO: Seperate model training/prediction from feat.extraction!'
    if GetTestFeatures==True:
        print('Starting to train model')
        #ORIG \\
##        model, lb_encoder = trainClassifier(trainingDir+'\\trainingSetFeatures.csv', False, 'forest', 0, False, False)
        model, lb_encoder = trainClassifier(filename=trainingDir+'/trainingSetFeatures.csv',normFlag= False,classifierType= 'forest',kbest= 0,alpha= False,optimalFlag= False)
        print('Created and trained model')


    if GetTestFeatures==True:
        print('Starting to extract features from testing set')
        #ORIG \\
##        featExt(testingDir, False, 'dir', trainingDir+'\\trainingSetNormParams.csv')
        featExt(testingDir, False, 'dir', trainingDir+'/trainingSetNormParams.csv')
        print('Extracted testing data features')
        #ORIG \\
##        df = pd.DataFrame.from_csv(testingDir+'\\testingSetFeatures.csv')
        df = pd.DataFrame.from_csv(testingDir+'/testingSetFeatures.csv')
        features = df.values
        print('Predicting labels')
        results = model.predict(features)

        'TODO: MAP this func, not if row by row! (SLOW!!)'
        ind = 0
        df['classname'] = 0
        for index, row in df.iterrows():
            df['classname'][index] = lb_encoder.inverse_transform(results[ind])
            ind+=1
        df = df['classname']
##        df.to_csv(testingDir+'\\testingSetResults.csv')
##        print('saved results to ' + testingDir+'\\testingSetResults.csv')
        df.to_csv(testingDir+'/testingSetResults.csv')
        print('saved results to ' + testingDir+'/testingSetResults.csv')
        writeClassifiedFastas(classType, testingDir, resultsDir, df)

if __name__ == '__main__' :
    pipeline()
