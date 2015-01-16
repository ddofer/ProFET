# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 18:24:07 2014

@author: Dan

Common tasks in processing data and models, likely to be called as part of the pipeline.
"""


from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier,AdaBoostClassifier,GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC,LinearSVC
from sklearn.naive_bayes import GaussianNB
from sklearn.grid_search import GridSearchCV
from sklearn.lda import LDA
from sklearn.decomposition import PCA
from operator import itemgetter
from Model_trainer import load_data
import numpy as np
from sklearn.preprocessing import StandardScaler

from sklearn.feature_selection import RFE, RFECV, FDR
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import RandomizedLogisticRegression

'''
https://pythonhosted.org/nolearn/_modules/nolearn/model.html#AveragingEstimator
Gets the majority vote for ensemble of classifiers
'''


'TODO. Partic featnames, and mrmr?'
def GetKFeatures(X,y,featNames, method='RFE',K=30, reduceMatrix = True):
    '''
    Gets best features using chosen method
    (K-best, RFE, RFECV,'L1' (RandomizedLogisticRegression),'Tree' (ExtraTreesClassifier), mrmr),
    then prints top K features' names (from featNames).
    If reduceMatrix =  True, then also returns X reduced to the K best features.

    Available methods' names are: 'RFE','RFECV','RandomizedLogisticRegression','K-best','ExtraTreesClassifier'..
    Note, that effectiveyl, Any scikit learn method could be used, if correctly imported..
    '''
    est = method()








def report(grid_scores, n_top=2) :
    '''
    Print out top models/parameters after a grid search for model params.
    '''
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores) :
        print("Model with rank: {0}".format(i + 1))
        print("Mean validation score: {0:.2f} (std: {1:.2f})".format(
            score.mean_validation_score, np.std(score.cv_validation_scores)))
        print("Parameters: {0}".format(score.parameters))
        print("")


def ModelParam_GridSearch(X_train, y_train, cv=4):
    '''
    Basic grid searchCV for multiple classifiers' perf & parameters.
    This is very limited and computationally expensive.
    Not guaranteed to reach even a local optima, but good to get a
    rough idea of parameters for the classifiers. (Does not address pre-processing)
    More classifiers can be added as desired, and parameters expanded.

    Later: Add options for RBM + Logit; PCA; LDA. See also
    http://scikit-learn-laboratory.readthedocs.org/en/latest/_modules/skll/learner.html

    TODO: Add parameters + put classifiers/"pipeline_#" in a list. (To allow checking only some params)
    '''

#    pipeline1 = Pipeline('clf', RandomForestClassifier() )
#
#    pipeline2 = Pipeline(
#    ('clf', KNeighborsClassifier()),)
    pipeline1 = RandomForestClassifier()
    pipeline2 = KNeighborsClassifier()
    pipeline3 = SVC()
    pipeline4 = GaussianNB()
    pipeline5 = AdaBoostClassifier()
    pipeline6 = SGDClassifier()
    pipeline7 = LogisticRegression()


    'RandomForestClassifier:'
    parameters1 = {
    'n_estimators': [85],
    'criterion': ['gini', 'entropy'],
    'max_features': ['auto',0.1,0.45],
    'max_depth': [11,None]
    }

    'KNeighborsClassifier:'
    parameters2 = {
    'n_neighbors': [3, 6, 9],
    'weights': ['uniform', 'distance']
    }

    'SVC:'
    parameters3 = {
    'C': [0.001,0.01, 0.1, 1.0,5],
    'kernel': ['rbf', 'poly','sigmoid'],
    'gamma': [0.01, 0.1, 1.0,5,20],
    'class_weight':['auto'], 'cache_size':[1500],
    }
    parameters4 = {}

    'AdaBoost:'
    parameters5 = {
    'n_estimators':[45,90], 'learning_rate':[0.01, 0.1,1]
    }

    'SGDClassifier:'
    parameters6 = {
     'alpha': [0.0001,0.00001,0.001],
    'penalty': ['l1','l2', 'elasticnet'],
    'n_iter': [60,160],
    'loss':['hinge','log', 'modified_huber'],'class_weight':['auto']}

    'LogisticRegression:'
    parameters7 = {
    'C': [0.001,.01, 0.1, 1.0,0.05,5,20],
    'penalty': ['l1','l2'],'class_weight':['auto']
    }

    'TODO: make this into a seperate method, with pars, pips passed to it as params'
    pars = [parameters1, parameters2, parameters3, parameters4,parameters5,parameters6,parameters7]
    pips = [pipeline1, pipeline2, pipeline3, pipeline4,pipeline5,pipeline6,pipeline7]

    print ("starting Gridsearch")
    for i in range(len(pars)):
        print(pips[i])
        gs = GridSearchCV(estimator=pips[i], param_grid=pars[i],
                          verbose=0, refit=False, n_jobs=2, cv=cv)
        gs = gs.fit(X_train, y_train)
        print ("Finished Gridsearch")
        print (gs.best_score_)
        report(gs.grid_scores_)

if __name__ == '__main__' :

    fileName = r'\trainingSetFeatures.csv'

    # filePath = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\CODE\Feature_Extract\test_seq\Chap'
    filePath = r'E:\Dropbox\Dropbox\BioInformatics Lab\AA_Information\FASTA_Sets\Benchmarks\LocTree3'

    # features, labels, lb_encoder,featureNames = load_data(filename, 'file')
    features, labels, lb_encoder,featureNames = load_data(filePath+fileName, 'file')

    X, y = features, labels
    print(X.shape,"X = samples, features")
    scale = StandardScaler(copy=False)
    X = scale.fit_transform(X)

    FD = SelectFdr( alpha=0.001)
    X = FD.fit_transform(X,y)
    print(X.shape,"X post FDR filt")

    pca = PCA(n_components=0.99)
    X_PCA = pca.fit_transform(X)
    print(X_PCA.shape,"X - PCA")
    print("X_PCA \n")
    ModelParam_GridSearch(X_PCA,y,cv=6)

    print("\n BASE X models: \n")
    ModelParam_GridSearch(X,y,cv=6)