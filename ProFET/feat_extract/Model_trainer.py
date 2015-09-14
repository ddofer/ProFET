__author__ = 'Michael'
'Modified from Michael doron - 7.8.2014'
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, \
    GradientBoostingClassifier
import numpy as np
from sklearn import metrics, cross_validation  #, linear_model, preprocessing)
from sklearn.svm import LinearSVC, SVC
from sklearn.metrics import precision_score, accuracy_score, recall_score
from sklearn.preprocessing import StandardScaler, LabelEncoder, \
    LabelBinarizer
from sklearn.lda import LDA
from sklearn import svm, preprocessing
from sklearn.linear_model import LogisticRegression, SGDClassifier

#from numba.decorators import jit
from time import time
from scipy.stats import randint as sp_randint
from sklearn.grid_search import RandomizedSearchCV,GridSearchCV
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import SelectPercentile, f_classif, SelectFwe, SelectFdr
from io import StringIO
from os.path import exists as file_exists

from Bio.SeqIO.FastaIO import SimpleFastaParser
from operator import itemgetter
from sklearn.pipeline import Pipeline

import pandas as pd

from IPython.core.debugger import Tracer #TO REMOVE!!!

def load_data(dataFrame="Feat_normalized.csv") :
    '''
    Load training data from csv file.  Load labels from it.
    Return matrix, training labels, encoder for labels.
    label_encoder uses transforming textual labels to integer and vice versa:
    label_encoder.inverse_transform(0) => 'Mammal_melanosome_0'
    label_encoder.transform('Mammal_melanosome_0') => 0
    '''
    label_exists = False
    #Tracer()() #TO REMOVE!!!
    if type(dataFrame) == type(''):
        if len(dataFrame) > 120 or not file_exists(dataFrame):
            #If it a string in format of csv and not a filename
            dataFrame = StringIO(dataFrame)
        #Load file
        df = pd.read_csv(dataFrame, delimiter='\t', header=0)
        try:
            df.set_index(keys = ['accession', 'classname'], inplace=True)
            #df = pd.read_csv(dataFrame, delimiter='\t', header=0, index_col=['accession', 'classname'])
            label_exists = True
        #When not labeled
        except KeyError:
            print('Features files does not contains labels')
            #df = pd.read_csv(dataFrame, delimiter='\t', header=0, index_col='accession')
            df.set_index(keys = 'accession', inplace=True)
    else:
        df = dataFrame

    features = df.values
    # M: creates numpy array
    feature_names=df.columns.values
    print("%s features" % (len(feature_names)))
    # print("feature_names: %s" %(feature_names))

    if label_exists:
        # create an object of scikit label encoder that transforms strings to ints
        label_encoder = LabelEncoder()
        # M: check if works with multiindex
        # M: take only label, not protein name

        s = df.index.get_level_values('classname')
        #print(s.value_counts)

        labels = label_encoder.fit_transform((df.index.get_level_values('classname').values))
        #print ("labels: %s %s" %(type(labels),labels))
        print("labels List: ",list(label_encoder.classes_))
        # M: creates numpy matrix (or array)
        return (features, labels, label_encoder, feature_names)
    accessions = df.index.get_level_values('accession')
    #To change the order!! (and in the calling functions)
    return features, accessions, feature_names

'TODO: Set params of model, by user selected/pretuned values (obtained in pipetakss via CV)'
def getClassifier(classifierType):
    if (classifierType == 'SGD'):
        model = SGDClassifier(penalty='elasticnet',class_weight='auto',n_jobs=-1,n_iter=150,l1_ratio =0.2)
    if (classifierType == 'LSVC'):
       model = LinearSVC(class_weight='auto')
    if (classifierType == 'forest'):
        model = RandomForestClassifier(n_jobs=-1, bootstrap=True, n_estimators=350,
                                        min_samples_leaf=1, min_samples_split =2,
                                        oob_score=False,max_features='auto',
                                        criterion='gini')
    if (classifierType == 'SVCrbf'):
       model = SVC(kernel="rbf", class_weight="auto", cache_size=1400, shrinking=True)
    if (classifierType == 'SVCpoly'):
        model = SVC(kernel="poly", class_weight="auto", cache_size=1400, shrinking=True)
    return model

def scale_data(X):
    '''Scale data by standard scaling
    Returns the scaled data and the scaler object
    '''
    normal_scaler = StandardScaler()
    X = normal_scaler.fit_transform(X)
    return X, normal_scaler

'TODO: Remove unneeded params here. MAke Callable. Add feat.sel by other methods (MI-mrmr ; L1-SVM, RFE..).'
'See port in PipeTasks.py - GetKFeatures '
def featureFitting(filename, X, y, featureNames,optimalFlag, kbest=20, alpha=0.05, model=None):
    '''
    Gets the K-best features (filtered by FDR, then select best ranked by t-test, more advanced options can be implemented).
    Save the data/matrix with the resulting/kept features to a new output file, "REDUCED_Feat.csv"
    Returns new features matrix, FD scaler, and K-select scaler
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
    return Reduced_df, FD, selectK

def trainClassifier(filename, classifierType, kbest, alpha, optimalFlag, normFlag=True):
    # extract the features, labels, lb encoder, and feature names
    features, labels, label_encoder, featureNames = load_data(filename)
    X, y = features, labels

    # change the names as ints back to strings
    class_names=label_encoder.inverse_transform(y)
    print('%s sequences, %s features each' % X.shape)

    if (normFlag == True):
        # 'Normalizing Unneeded if features already normalized. (scaling may still be needed)'
        X, scaler = scale_data(X)
        print("Features Data scaled")

    # create classifier model object
    model = getClassifier(classifierType)

    # if needed, fit the features
    if ((kbest != 0) or (alpha == True) or (optimalFlag == True)):
        df, fd_scaler, kselect_scaler = featureFitting(filename=filename, X=X, y=y, kbest=kbest, featureNames=featureNames, alpha=alpha, optimalFlag=optimalFlag,model=model)
        features, labels, label_encoder,featureNames = load_data(df)
        X, y = features, labels
    model.fit(X, y)
    return model, label_encoder, scaler, featureNames



if __name__ == '__main__' :
    trainClassifier()
