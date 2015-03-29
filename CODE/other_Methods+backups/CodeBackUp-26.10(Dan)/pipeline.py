#!/sw/bin/python3.3
#! E:\Python33\python

__author__ = 'Michael'
'Modified from Michael doron - 7.8.2014'
from FeatureGen import featExt
from Model_trainer import trainClassifier

import time
import pandas as pd
import sklearn
import numpy
from sklearn.preprocessing import MinMaxScaler, StandardScaler, LabelEncoder

import click

# trainD=r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\CODE\Feature_Extract\test_seq\Chap'
# testD=r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\FASTA_Sets\HSP33_Chap\Unknown_Tests.2'

# trainD=r'D:\SkyDrive\Dropbox\BioInformatics Lab\AA_Information\CODE\Feature_Extract\test_seq\Chap'
# trainD=r'D:\SkyDrive\Dropbox\BioInformatics Lab\AA_Information\CODE\Feature_Extract\test_seq\TestFilt'
trainD=r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\FASTA_Sets\HSP33_Chap\V4_October'

testD=r'D:\SkyDrive\Dropbox\BioInformatics Lab\AA_Information\FASTA_Sets\HSP33_Chap\Unknown_Tests.2'

@click.command()
@click.option('--trainingSetDir','-r','trainingDir',default=trainD,
              help='The path to the training set fasta files', type = str)
@click.option('--testingSetDir','-s','testingDir',default=testD,
              help='The path to the testing set fasta files', type = str)
@click.option('--trainFeatures','-rf','GetTrainingFeatures',default=True,help='Whether to extract the training set features', type = bool)
@click.option('--testFeatures','-sf','GetTestFeatures',default=True,help='Whether to get the testing set features', type = bool)
@click.option('--classType','-ct','classType',default='file',help="Defines the classname of each protein, by \'dir\', \'file\', or \'id\'.", type = str)
# @click.option('--classType','-ct','classType',default='id',help="Defines the classname of each protein, by \'dir\', \'file\', or \'id\'.", type = str)

def pipeline(trainingDir,testingDir,GetTrainingFeatures,GetTestFeatures, classType):
    # change here to the training data folder
    # trainingDir = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\CODE\Feature_Extract\test_seq\Chap'
    # trainingDir = str(trainingDir)

    # change here to the testing data folder
    # testingDir = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\FASTA_Sets\HSP33_Chap\Unknown_Tests'
    # testingDir = str(testingDir)

    if GetTrainingFeatures==True:
        print('Starting to extract features from training set')
        'Temporary measure: If features extracted and saved, disable following line to avoid re-extracting trainign features'
        t0 = time.clock()
        featExt(directory=trainingDir, trainingSetFlag=True,
                classType=classType, normFlag=True, normParams='.')
        t1 = time.clock()
        time1 = t1-t0
        print('Extracted training data features')
	    #

    print('Starting to train model')
    t1 = time.clock()
    model, lb_encoder = trainClassifier(trainingDir+'\\trainingSetFeatures.csv', False, 'forest', 0, False, False)
    t1 = time.clock()
    time2 = t1-t0
    print('Created and trained model')


    if GetTestFeatures==True:
        print('Starting to extract features from testing set')
        t0 = time.clock()
        featExt(directory=testingDir, trainingSetFlag=False, classType='dir', normFlag=True,normParams= trainingDir+'\\trainingSetNormParams.csv')
        t1 = time.clock()
        time3 = t1-t0
        print('Extracted testing data features')
        df = pd.DataFrame.from_csv(testingDir+'\\testingSetFeatures.csv')
        features = df.values
        print('Predicting labels')
        t0 = time.clock()
        results = model.predict(features)
        t1 = time.clock()
        time4 = t1-t0
        ind = 0

        df['classname'] = 0
        for index, row in df.iterrows():
            df['classname'][index] = lb_encoder.inverse_transform(results[ind])
            ind+=1
        df = df['classname']
        df.to_csv(testingDir+'\\testingSetResults.csv')
        print('saved results to ' + testingDir+'\\testingSetResults.csv')

if __name__ == '__main__' :
    # cType = 'file'
    #Not Working - Dan #
#    pipeline(trainingDir=testD,testingDir=testD,GetTrainingFeatures=True,GetTestFeatures=True,classType=cType)

    pipeline(GetTestFeatures=False)

    # featExt
