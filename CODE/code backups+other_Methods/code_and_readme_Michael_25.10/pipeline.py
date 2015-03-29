#!/sw/bin/python3.3
#! E:\Python33\python

__author__ = 'Michael'
'Modified from Michael doron - 7.8.2014'


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
@click.option('--testFeatures','-sf','GetTestFeatures',default=True,help='Whether to get the testing set features', type = bool)
@click.option('--classType','-ct','classType',default='id',help="Defines the classname of each protein, by \'dir\', \'file\', or \'id\'.", type = str)
def pipeline(trainingDir,testingDir,resultsDir, GetTrainingFeatures,GetTestFeatures, classType):
    # change here to the training data folder
    # trainingDir = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\CODE\Feature_Extract\test_seq\Chap'
    # trainingDir = str(trainingDir)

    # change here to the testing data folder
    # testingDir = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\FASTA_Sets\HSP33_Chap\Unknown_Tests'
    # testingDir = str(testingDir)

    # change here to the type of classification (by 'file' or by 'dir')

    if GetTrainingFeatures==True:
        print('Starting to extract features from training set')
        'Temporary measure: If features extracted and saved, disable following line to avoid re-extracting trainign features'
        featExt(directory=trainingDir, trainingSetFlag=True,
                classType=classType, normParams='.')
        print('Extracted training data features')
	    #

    print('Starting to train model')
    model, lb_encoder = trainClassifier(trainingDir+'\\trainingSetFeatures.csv', False, 'forest', 0, False, False)
    print('Created and trained model')


    if GetTestFeatures==True:
        print('Starting to extract features from testing set')
        featExt(testingDir, False, 'dir', trainingDir+'\\trainingSetNormParams.csv')
        print('Extracted testing data features')
        df = pd.DataFrame.from_csv(testingDir+'\\testingSetFeatures.csv')
        features = df.values
        print('Predicting labels')
        results = model.predict(features)
        ind = 0

        df['classname'] = 0
        for index, row in df.iterrows():
            df['classname'][index] = lb_encoder.inverse_transform(results[ind])
            ind+=1
        df = df['classname']
        df.to_csv(testingDir+'\\testingSetResults.csv')
        print('saved results to ' + testingDir+'\\testingSetResults.csv')
        writeClassifiedFastas(classType, testingDir, resultsDir, df)

if __name__ == '__main__' :
    pipeline()
