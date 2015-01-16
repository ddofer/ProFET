#!/sw/bin/python3.3
#! E:\Python33\python

__author__ = 'Michael'
'Modified from Michael doron - 7.8.2014'

from Feature_Extract.FeatureGen import featExt
from Machine_Learn.Model_trainer import trainClassifier

# from FeatureGen import featExt
# from Model_trainer import trainClassifier

import pandas as pd
import sklearn
import numpy
from sklearn.preprocessing import MinMaxScaler, StandardScaler, LabelEncoder

def pipeline():
    # change here to the training data folder
    # trainingDir = 'C:\\Users\\Michael\\Desktop\\DanCode\\Feature_Extract\\test_seq\\example\\NP+SP+Negs_BigBalanced'
    trainingDir = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\CODE\Feature_Extract\test_seq\Chap'

    # change here to the testing data folder
    # testingDir = 'C:\\Users\\Michael\\Desktop\\DanCode\\Feature_Extract\\test_seq\\example\\NP+SP+Negs_smallBalanced'
    testingDir = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\FASTA_Sets\HSP33_Chap\Unknown_Tests'

    # change here to the type of classification (by 'file' or by 'dir')
    classificationType = 'file'
    print('Starting to extract features from training set')
    featExt(trainingDir, True, classificationType, True, '.')
    print('Extracted training data features')
    print('Starting to train model')
    model, lb_encoder = trainClassifier(trainingDir+'\\trainingSetFeatures.csv', False, 'forest', 0, False, False)
    print('created and trained model')
    print('starting to extract features from testing set')
    featExt(testingDir, False, 'dir', True, trainingDir+'\\trainingSetNormParams.csv')
    print('extracted testing data features')
    df = pd.DataFrame.from_csv(testingDir+'\\testingSetFeatures.csv')
    features = df.values
    print('predicting the labels')
    results = model.predict(features)
    ind = 0
    df['classname'] = 0
    for index, row in df.iterrows():
        df['classname'][index] = lb_encoder.inverse_transform(results[ind])
        ind+=1
    df = df['classname']
    df.to_csv(testingDir+'\\testingSetResults.csv')
    print('saved results to ' + testingDir+'\\testingSetResults.csv')

if __name__ == '__main__' :
    pipeline()
