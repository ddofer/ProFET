"""
Given a CSV file with features - get it's performance/results!
This is meant to be used after we've settled on hyperparameters, features, etc'.
It's to get "scores" for quoting in the article. (barring the test set)
"""

#http://pythonhosted.org/mlxtend/#ensembleclassifier
from ensemble import EnsembleClassifier

import os
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, \
    GradientBoostingClassifier,AdaBoostClassifier
import numpy as np
from sklearn import metrics, cross_validation  #, linear_model, preprocessing)
from sklearn.svm import LinearSVC, SVC
from sklearn.metrics import precision_score, accuracy_score, recall_score
from sklearn import svm, preprocessing
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression,RandomizedLogisticRegression,LogisticRegressionCV
from sklearn import metrics
from scipy.stats import randint as sp_randint
from sklearn.grid_search import RandomizedSearchCV,GridSearchCV
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif, chi2, \
    SelectFwe, SelectFdr
from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import StratifiedKFold,cross_val_score,StratifiedShuffleSplit,cross_val_predict
from sklearn.naive_bayes import MultinomialNB,GaussianNB
from sklearn.utils import shuffle
from sklearn.calibration import CalibratedClassifierCV
import pandas as pd

from api.ml import *
from api import config
from api.classifier_training import get_balanced_weights
from sklearn.ensemble import BaggingClassifier
from sklearn.neighbors import KNeighborsClassifier
from api.classifier_training import get_training_data

# Use a constant seed
np.random.seed(1274)
SEED = 14  # always use a seed for randomized procedures

"""
Feature selection should be done within a CV loop, not externally.
We test if Feature selection is good with gridsearch+pipeline.

"""

"""
We may use an ensemble of models (implemented in mlextend - downloaded locally to ensemble.py):
http://nbviewer.ipython.org/github/rasbt/mlxtend/blob/master/docs/examples/sklearn_ensemble_ensembleclassifier.ipynb
"""

def get_scores(scores,y,label):
    '''
    Returns a dictionary of metrics for a given classification of the data (given by Cross_val_predict).
    scores: list
        Classifier predictions on data
    y: list
        True Class labels
    label: string
        Name of the classifier used
    '''

    roc_auc = metrics.roc_auc_score(y, scores,average=None)
    print("roc_auc (No-Av): %0.4f " % (roc_auc))
    roc_auc = metrics.roc_auc_score(y, scores,average='weighted')
    print("roc_auc (weighted-Av): %0.4f " % (roc_auc))
    f1_pos = metrics.f1_score(y, scores,average='binary')
    print("POS f1: %0.4f  " % (f1_pos))
    av_PR = metrics.average_precision_score(y, scores) # corresponds to the area under the precision-recall curve
    print("Av_precision (Prec-Recall AUC): %0.3f " % (av_PR))
    accuracy = metrics.accuracy_score(y, scores)
    print("Accuracy: %0.3f " % (accuracy))
    precision,recall,fscore,support = metrics.precision_recall_fscore_support(y, scores,average='binary')
    print("Precision: %0.3f " % (precision))
    print("Recall: %0.3f " % (recall))
    # print("fscore(fBeta): %0.4f  [%s]" % (fscore, label))
    mcc = metrics.matthews_corrcoef(y, scores)
    print("MCC: %0.3f " % (mcc))

    results_dict = {'roc_auc(macro)':roc_auc,'f1_Pos':f1_pos,'accuracy':accuracy,
    'precision':precision,'recall':recall,
    # 'fscore-fBeta':fscore,
    'average Precision':av_PR,'mcc':mcc
    }
    results_dict = {k:round(float(v),4) for k, v in results_dict.items()}


    return results_dict



# E:\Dropbox\Dropbox\Protein Cleavage Prediction\data\NeuroPred\V4\features-11_8_KR.csv
def testEnsemble(target_file_loc = r'D:\Dropbox\Protein Cleavage Prediction\data\NeuroPred\V4\best_windows_pos\features10_8.csv',outputFileName="MetricsResults"):
    """
    http://nbviewer.ipython.org/github/rasbt/mlxtend/blob/master/docs/examples/sklearn_ensemble_ensembleclassifier.ipynb#Additional-Note-About-the-EnsembleClassifier-Implementation:-Class-Labels-vs.-Probabilities
    http://sebastianraschka.com/Articles/2014_ensemble_classifier.html#EnsembleClassifier---Tuning-Weights

    Blend:
    https://github.com/log0/vertebral/blob/master/stacked_generalization.py
    """

    np.random.seed(123)

    SILLY_NUMBER = 70 #USed as a magic number; when debugging for speed

    results = {} #List of classifier results-dicts.

    clf1 = LogisticRegressionCV(Cs=22,class_weight='auto')
    clf2 = RandomForestClassifier(n_estimators=int(SILLY_NUMBER*1.5), max_features=SILLY_NUMBER,bootstrap=False,class_weight='auto',n_jobs=2,criterion='entropy',random_state=123)
    clf3 = SVC(C=2.3, kernel= 'rbf', gamma= 0.0, cache_size= 1000, class_weight= 'auto',probability=True)
    # C=3.798 , 5.79
    clf4 = GradientBoostingClassifier(n_estimators=SILLY_NUMBER, max_depth=12,min_samples_leaf=2)
    clf5 = BaggingClassifier(KNeighborsClassifier(), max_samples=0.6, max_features=0.4)
    # clf6 = Pipeline([('scale',MinMaxScaler(copy = False)),('MultinomialNB',MultinomialNB())])
    clf7 = KNeighborsClassifier(n_neighbors=4, weights='distance')
    clf8 = SVC(C=20,class_weight= 'auto',probability=True)

    cclf1,cclf2,cclf3,cclf8 = CalibratedClassifierCV(clf1),CalibratedClassifierCV(clf2,cv=4),CalibratedClassifierCV(clf3,cv=4),CalibratedClassifierCV(clf8,cv=4)
    clfs_calibrated = [cclf1,cclf2,cclf3]
    cclf5,cclf7 = CalibratedClassifierCV(clf5),CalibratedClassifierCV(clf7)


    X,y,KM_pred = get_training_data(target_file_loc, drop_duplicates = True, select_features = True, scale = True,get_KM = True)

    # eclf = EnsembleClassifier(clfs=[clf1, clf2, clf3,clf4, clf5], voting='soft')#, weights=y_weights)
    eclf = EnsembleClassifier(clfs=[clf1, cclf2, cclf3,cclf8,cclf7], voting='soft', weights=[2.5,1,2,1])

    eclf2 = EnsembleClassifier(clfs=[clf1,cclf2,clf3,clf4,clf8], voting='hard')#,weights=[2,2,2.5,1,1])


    all_clfs_calibrated = [cclf1,cclf2,cclf3,cclf5,cclf7,cclf8,eclf]
    # eclf_all = EnsembleClassifier(clfs=all_clfs_calibrated, voting='hard')

    classifiers_and_names = zip([
        # clf1, cclf2,
         clf3,
     # cclf5,cclf7,cclf8,
     eclf,
     eclf2,
     # eclf_all
     ],
     [
     # 'Logistic Regression','Random Forest',
     'SVM-RBF',
     # 'BaggingClassifier-KNN','KNeighbors',  'linearSVC',
      'Ensemble-SoftWeighted','Ensemble-Hard',
      # 'Ensemble-Ensemble-all'
      ])

    # classifiers_and_names = zip([eclf,eclf2],['Ensemble-SoftWeighted','Ensemble-Hard'])

    print("X Shape:",X.shape)
    print('# Positives: %i' % (sum(y)))

    for clf, label in classifiers_and_names:
        print('')
        print(label)
        print('')
        scores = cross_val_predict(clf, X, y, cv=8,n_jobs=-1)
        results[label] = get_scores(scores,y,label)



    print('Predicted according to Known Motif model: %i' %(sum(KM_pred)))
    cm = metrics.classification_report(y, KM_pred)
    print(cm)
    cm = metrics.confusion_matrix(y, KM_pred)
    print(cm)
    results['KnownMotif'] = get_scores(KM_pred,y,label='KnownMotif')

    res_df = pd.DataFrame(results)
    res_df.to_csv(outputFileName+"tsv", sep='\t')
    res_df.to_csv(outputFileName+"csv", sep=';')



def blend(target_file_loc = r'E:\Dropbox\Dropbox\Protein Cleavage Prediction\data\NeuroPred\V4\features-11_8_KR.csv'):
    '''
    https://github.com/log0/vertebral/blob/master/stacked_generalization.py#L162
    '''
    X,Y = get_training_data(target_file_loc, drop_duplicates = True, select_features = True, scale = True)
    # We need to transform the string output to numeric
    # label_encoder = LabelEncoder()
    # label_encoder.fit(Y)
    # Y = label_encoder.transform(Y)

    # The DEV SET will be used for all training and validation purposes
    # The TEST SET will never be used for training, it is the unseen set.
    dev_cutoff = len(Y) * 4/5
    X_dev = X[:dev_cutoff]
    Y_dev = Y[:dev_cutoff]
    X_test = X[dev_cutoff:]
    Y_test = Y[dev_cutoff:]

    n_trees = 30
    n_folds = 4

    # Our level 0 classifiers
    clfs = [
        ExtraTreesClassifier(n_estimators = n_trees * 2, criterion = 'gini'),
        LogisticRegressionCV(Cs=25,class_weight='auto'),
    RandomForestClassifier(n_estimators=150, max_features=100,bootstrap=False,class_weight='auto',n_jobs=-2,criterion='entropy',random_state=123),
    SVC(C=3.798, kernel= 'rbf', gamma= 0.0, cache_size= 1000, class_weight= 'auto',probability=True),
    GradientBoostingClassifier(n_estimators=110, max_depth=9,min_samples_leaf=2),
    BaggingClassifier(KNeighborsClassifier(), max_samples=0.6, max_features=0.5),
    # Pipeline([('scale',MinMaxScaler(copy = False)),('MultinomialNB',MultinomialNB())]),
     KNeighborsClassifier(n_neighbors=5, weights='distance')
    ]

    # Ready for cross validation
    skf = list(StratifiedKFold(Y_dev, n_folds))

    # Pre-allocate the data
    blend_train = np.zeros((X_dev.shape[0], len(clfs))) # Number of training data x Number of classifiers
    blend_test = np.zeros((X_test.shape[0], len(clfs))) # Number of testing data x Number of classifiers

    print ('X_test.shape = %s' % (str(X_test.shape)))
    print ('blend_train.shape = %s' % (str(blend_train.shape)))
    print ('blend_test.shape = %s' % (str(blend_test.shape)))

    # For each classifier, we train the number of fold times (=len(skf))
    for j, clf in enumerate(clfs):
        print ('Training classifier [%s]' % (j))
        blend_test_j = np.zeros((X_test.shape[0], len(skf))) # Number of testing data x Number of folds , we will take the mean of the predictions later
        for i, (train_index, cv_index) in enumerate(skf):
            print ('Fold [%s]' % (i))

            # This is the training and validation set
            X_train = X_dev[train_index]
            Y_train = Y_dev[train_index]
            X_cv = X_dev[cv_index]
            Y_cv = Y_dev[cv_index]

            clf.fit(X_train, Y_train)

            # This output will be the basis for our blended classifier to train against,
            # which is also the output of our classifiers
            ## Orig
            # blend_train[cv_index, j] = clf.predict(X_cv)
            # blend_test_j[:, i] = clf.predict(X_test)
            blend_train[cv_index, j] = clf.predict_proba(X_cv)[:,1]
            blend_test_j[:, i] = clf.predict_proba(X_test)[:,1]
        # Take the mean of the predictions of the cross validation set
        blend_test[:, j] = blend_test_j.mean(1)

    print ('Y_dev.shape = %s' % (Y_dev.shape))

    # Start blending!
    # bclf = LogisticRegressionCV()
    bclf = AdaBoostClassifier(n_estimators=60)
    bclf.fit(blend_train, Y_dev)

    # Predict now
    Y_test_predict = bclf.predict(blend_test)
    score = metrics.accuracy_score(Y_test, Y_test_predict)
    print ('Accuracy = %s' % (score))
    score = metrics.f1_score(Y_test, Y_test_predict)
    print ('f1 = %s' % (score))
    score = metrics.roc_auc_score(Y_test, Y_test_predict)
    print ('roc_auc = %s' % (score))
    return score


##############################################################################

if __name__ == '__main__' :
    TARGET_DIR = r'D:\Dropbox\Protein Cleavage Prediction\data\NeuroPred\V4\altFeats'
    #
    os.chdir(TARGET_DIR)
    # target_files = ['features10_9_KR_KM.csv', 'features10_9_KR.csv', 'features11_8_KR.csv', 'features11_8_KR_KM.csv']
    target_files=['features_redAA.csv']

    for f in target_files:
        print("File:", f)
        fname = "Metrics-"+(f.replace("csv","").replace("features",""))
        testEnsemble(target_file_loc=f, outputFileName=str(fname))
        print

    # blend()
