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
from sklearn.feature_selection import SelectPercentile, f_classif, SelectFwe, SelectFdr

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

    s=df.index.get_level_values('classname')
    #print(s.value_counts)

    labels = lb.fit_transform((df.index.get_level_values('classname').values))
    #print ("labels: %s %s" %(type(labels),labels))
    print("labels List: ",list(lb.classes_))
    # M: creates numpy matrix (or array)
    features = df.values
    # M: creates numpy array
    feature_names=df.columns.values
    print("%s features" % (len(feature_names)))
    # print("feature_names: %s" %(feature_names))
    return (features, labels, lb,feature_names)

'TODO: Set params of model, by user selected/pretuned values (obtained in pipetakss via CV)'
def getClassifier(classifierType):
    if (classifierType == 'SGD'):
        model = SGDClassifier(penalty='elasticnet',class_weight='auto',n_jobs=-1,n_iter=75,l1_ratio =0.2)
    if (classifierType == 'LSVC'):
       model = LinearSVC(class_weight='auto')
    if (classifierType == 'forest'):
        model = RandomForestClassifier(n_jobs=-1, bootstrap=True, n_estimators=500,
                                        min_samples_leaf=1, min_samples_split =2,
                                        oob_score=True,max_features='auto',
                                        criterion='gini')
    if (classifierType == 'SVCrbf'):
       model = SVC(kernel="rbf", class_weight="auto", cache_size=1400, shrinking=True)
    if (classifierType == 'SVCpoly'):
        model = SVC(kernel="poly", class_weight="auto", cache_size=1400, shrinking=True)
    return model


def scale_data(X) :
    scale_scaler = MinMaxScaler()
    normal_scaler = StandardScaler()
    X = normal_scaler.fit_transform(X)
    X = scale_scaler.fit_transform(X)
    return X

'TODO: Remove unneeded params here. MAke Callable. Add feat.sel by other methods (MI-mrmr ; L1-SVM, RFE..).'
'See port in PipeTasks.py - GetKFeatures '
def featureFitting( filename, X, y, featureNames,optimalFlag, kbest=20, alpha=0.05,model=None):
    '''
    Gets the K-best features (filtered by FDR, then select best ranked by t-test , more advanced options can be implemented).
    Save the data/matrix with the resulting/kept features to a new output file, "REDUCED_Feat.csv"
    '''
    a=alpha
    FD = SelectFdr(alpha=a)
    X = FD.fit_transform(X,y)

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
        df = featureFitting( filename=filename, X=X, y=y, kbest=kbest, featureNames=featureNames, alpha=alpha, optimalFlag=optimalFlag,model=model)
        features, labels, lb_encoder,featureNames = load_data(df, 'df')
        X, y = features, labels
    model.fit(X, y)
    return model, lb_encoder



if __name__ == '__main__' :
    trainClassifier()
