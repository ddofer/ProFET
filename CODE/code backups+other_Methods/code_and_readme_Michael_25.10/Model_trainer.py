__author__ = 'Michael'
'Modified from Michael doron - 7.8.2014'
import click
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, \
    GradientBoostingClassifier
import numpy as np
from sklearn import metrics, cross_validation  #, linear_model, preprocessing)
from sklearn.svm import LinearSVC, SVC
from sklearn.metrics import precision_score, accuracy_score, recall_score
from sklearn.preprocessing import MinMaxScaler, StandardScaler, LabelEncoder, \
    LabelBinarizer
from sklearn.lda import LDA
from sklearn import svm, preprocessing
from sklearn.linear_model import LogisticRegression, SGDClassifier

from numba.decorators import jit
from time import time
from scipy.stats import randint as sp_randint
from sklearn.grid_search import RandomizedSearchCV,GridSearchCV
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import SelectPercentile, f_classif, chi2, \
    SelectFwe, SelectFdr

from Bio.SeqIO.FastaIO import SimpleFastaParser
from operator import itemgetter
from sklearn.pipeline import Pipeline

import pandas as pd



def load_data(dataFrame="Feat_normalized.csv", dfType='file') :
    '''
    Load training data from csv file.  Load labels from it.
    Return matrix, training labels, encoder for labels.
    '''
    if dfType == 'file':
        df = pd.read_csv(dataFrame, index_col=[0,1]) # is index column 0 in multiindex as well?
    else:
        df = dataFrame
    # create an object of scikit label encoder that transforms strings to ints
    lb = LabelEncoder()
    # M: check if works with multiindex
    # M: take only label, not protein name
    print (df.index.get_level_values('classname').values)
    labels = lb.fit_transform((df.index.get_level_values('classname').values))
    print ("labels: %s %s" %(type(labels),labels))
    # M: creates numpy matrix (or array)
    features = df.values
    # M: creates numpy array
    feature_names=df.columns.values
    print("%s features: " % (len(feature_names)))
    # classes = label_encoder.transform(np.asarray(df['labels']))
    print("labels: ")
    print(labels)
    # print("feature_names: %s" %(feature_names))
    return (features, labels, lb,feature_names)


def getClassifier(classifierType):
    if (classifierType == 'SGD'):
        model = SGDClassifier(penalty='elasticnet',class_weight='auto',n_jobs=-1,n_iter=35,l1_ratio =0.2)
    if (classifierType == 'LSVC'):
       model = LinearSVC(class_weight='auto')
    if (classifierType == 'forest'):
        model = RandomForestClassifier(n_jobs=-1, bootstrap=True, n_estimators=380,
                                        min_samples_leaf=2, min_samples_split =2,
                                        oob_score=True,max_features='auto',
                                        criterion='gini', max_depth=12)
    if (classifierType == 'SVCrbf'):
       model = SVC(kernel="rbf", class_weight="auto", cache_size=1200, shrinking=True)
    if (classifierType == 'SVCpoly'):
        model = SVC(kernel="poly", cache_size=1200, shrinking=True)
    return model


def scale_data(X) :
    scale_scaler = MinMaxScaler()
    normal_scaler = StandardScaler()
    X = normal_scaler.fit_transform(X)
    X = scale_scaler.fit_transform(X)
    return X


def featureFitting(model, filename, X, y, featureNames, kbest, alpha, optimalFlag):
    selectK = SelectKBest(k=kbest)
    selectK.fit(X,y)
    selectK_mask=selectK.get_support()
    K_featnames = featureNames[selectK_mask]
    print("K_featnames: %s" %(K_featnames))
    Reduced_df = pd.read_csv(filename, index_col=0)
    Reduced_df = Reduced_df[Reduced_df.columns[selectK_mask]]
    Reduced_df.to_csv('REDUCED_Feat.csv')
    return Reduced_df
#
# @click.command()
# @click.option('--training_set_file','-tsf','filename',default='.',help='The path to the trainings set dataframe csv')
# @click.option('--normalize','-n','normFlag',default=False,help='A flag indicating whether to normalize the data', type=bool)
# @click.option('--classifier','-c','classifierType',default='forest',help='The type of the classifier')
# @click.option('--kbest','-k','kbest',default=0,help='The number of features to use, chosen by k-best. enter 0 to not use kbest')
# @click.option('--alpha','-a','alpha',default=False,help='A flag indicating whether to use alpha feature fitting', type=bool)
# @click.option('--optimal','-o','optimalFlag',default=False,help='A flag indicating whether to use optimal feature transform', type=bool)
def trainClassifier(filename, normFlag, classifierType, kbest, alpha, optimalFlag):
    # extract the features, labels, lb encoder, and feature names
    features, labels, lb_encoder,featureNames = load_data(filename, 'file')
    X, y = features, labels

    # change the names as ints back to strings
    class_names=lb_encoder.inverse_transform(y)
    print("Data and labels imported")
    print(X.shape)

    if (normFlag == True):
        # 'Normalizing Unneeded if features already normalized. (scaling may still be needed)'
        X = scale_data(X)
        print("Features Data scaled")

    # create classifier model object
    model = getClassifier(classifierType)

    # if needed, fit the features
    if ((kbest != 0) or (alpha == True) or (optimalFlag == True)):
        df = featureFitting(model, filename, X, y, kbest, featureNames, alpha, optimalFlag)
        features, labels, lb_encoder,featureNames = load_data(df, 'df')
        X, y = features, labels
    model.fit(X, y)
    return model, lb_encoder



if __name__ == '__main__' :
    trainClassifier()
