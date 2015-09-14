# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 18:24:07 2014

@author: Dan

Common tasks in processing data and models, likely to be called as part of the pipeline.
"""

from sys import argv
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier,AdaBoostClassifier,GradientBoostingClassifier, ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier, RandomizedLogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC,LinearSVC,NuSVC
from sklearn.naive_bayes import GaussianNB
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import StratifiedKFold,cross_val_score,StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE, RFECV, SelectFdr,f_classif,SelectFwe,SelectPercentile,SelectKBest
from sklearn.linear_model import RandomizedLogisticRegression
from sklearn.lda import LDA
from sklearn.decomposition import PCA,FastICA #,TruncatedSVD
from operator import itemgetter
from Model_trainer import load_data
from sklearn.metrics import confusion_matrix, precision_recall_fscore_support,matthews_corrcoef, classification_report
#http://tokestermw.github.io/posts/imbalanced-datasets-random-forests/

# from sklearn.preprocessing import balance_weights

import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline

'''
https://pythonhosted.org/nolearn/_modules/nolearn/model.html#AveragingEstimator
Gets the majority vote for ensemble of classifiers

Good enhanced confusion matrix + weights for RF:
http://stackoverflow.com/questions/24123498/recursive-feature-elimination-on-random-forest-using-scikit-learn
'''

'Use SKLL? https://skll.readthedocs.org/en/latest/run_experiment.html#quick-example'


def balance_weights(y):
    """
    https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/preprocessing/_weights.py
    http://stackoverflow.com/questions/20082674/unbalanced-classification-using-randomforestclassifier-in-sklearn

    Compute sample weights such that the class distribution of y becomes
       balanced.
    Parameters
    ----------
    y : array-like
        Labels for the samples.
    Returns
    -------
    weights : array-like
        The sample weights.
    """
    y = np.asarray(y)
    y = np.searchsorted(np.unique(y), y)
    bins = np.bincount(y)

    weights = 1. / bins.take(y)
    weights *= bins.min()

    return weights

'TODO: Implement.  (More plots and test figures in source)'
def Feature_Importance_plot (est,names):
    'http://nbviewer.ipython.org/github/pprett/pydata-gbrt-tutorial/blob/master/gbrt-tutorial.ipynb'
    fx_imp = pd.Series(est.feature_importances_, index=names)
    fx_imp /= fx_imp.max()  # normalize
    fx_imp.sort()
    fx_imp.plot(kind='barh', figsize=FIGSIZE)

def createImportancePlot(splt,desc,importances,caption):
    '''
    http://nbviewer.ipython.org/gist/greglandrum/4316460
    '''
    import numpy.numarray as na
    labels = []
    weights = []
    threshold = sort([abs(w) for w in importances])[-11]
    for d in zip(desc,importances):
        if abs(d[1]) > threshold:
            labels.append(d[0])
            weights.append(d[1])
    xlocations = na.array(range(len(labels)))+0.5
    width = 0.8
    splt.bar(xlocations, weights, width=width)
    splt.set_xticks([r+1 for r in range(len(labels))])
    splt.set_xticklabels(labels)
    splt.set_xlim(0, xlocations[-1]+width*2)
    splt.set_title(caption)
    splt.get_xaxis().tick_bottom()
    splt.get_yaxis().tick_left()

def get_enhanced_confusion_matrix(actuals, predictions, labels):
    """"
    Enhances confusion_matrix by adding sensivity and specificity metrics
    http://stackoverflow.com/questions/24123498/recursive-feature-elimination-on-random-forest-using-scikit-learn
    """
    cm = confusion_matrix(actuals, predictions, labels = labels)
    sensitivity = float(cm[1][1]) / float(cm[1][0]+cm[1][1])
    specificity = float(cm[0][0]) / float(cm[0][0]+cm[0][1])
    weightedAccuracy = (sensitivity * 0.9) + (specificity * 0.1)
    return cm, sensitivity, specificity, weightedAccuracy

def Get_yPred (X,y,clf_class,n_folds=10, pred_proba=False) : #,**kwargs):
    '''
    Return "Full" Y_predictions from a given c;assifier (not just from one split): (From def run_cv)
    http://blog.yhathq.com/posts/predicting-customer-churn-with-sklearn.html

    Could also be done with stratified shuffle split (+Append output) ?
    http://scikit-learn.org/stable/modules/generated/sklearn.cross_validation.StratifiedShuffleSplit.html
    '''
    # Construct a kfolds object
    # kf = StratifiedKFold(len(y),n_folds,shuffle=True) #shuffle?
    kf = StratifiedKFold(y,n_folds,shuffle=True) #shuffle?
    y_pred = y.copy()

    # Iterate through folds
    for train_index, test_index in kf:
        X_train, X_test = X[train_index], X[test_index]
        y_train = y[train_index]
        # sample_weight=balance_weights(y_train)

        # Initialize a classifier with key word arguments
        clf = clf_class #(**kwargs)
        #sample_weight weighting not working here.. ?  TODO
        clf.fit(X_train,y_train) #,sample_weight) #
        if pred_proba == True:
            y_pred[test_index] = clf.predict_proba(X_test)
        else:
            y_pred[test_index] = clf.predict(X_test)
    return y_pred

'TODO. Partic featnames, and mrmr?'
'Code near duplicate of that in Model_train.py - def featureFitting'
'Currently works on the file, rather than dataframe/inmemory'
def GetKFeatures(filename, method='RFE',kbest=30,alpha=0.01, reduceMatrix = True):
    '''
    Gets best features using chosen method
    (K-best, RFE, RFECV,'L1' (RandomizedLogisticRegression),'Tree' (ExtraTreesClassifier), mrmr),
    then prints top K features' names (from featNames).
    If reduceMatrix =  True, then also returns X reduced to the K best features.

    Available methods' names are: 'RFE','RFECV','RandomizedLogisticRegression','K-best','ExtraTreesClassifier'..
    Note, that effectiveyl, Any scikit learn method could be used, if correctly imported..
    '''
    #est = method()
    '''
    Gets the K-best features (filtered by FDR, then select best ranked by t-test , more advanced options can be implemented).
    Save the data/matrix with the resulting/kept features to a new output file, "REDUCED_Feat.csv"
    '''
    features, labels, lb_encoder,featureNames = load_data(filename)
    X, y = features, labels

    # change the names as ints back to strings
    class_names=lb_encoder.inverse_transform(y)
    print("Data and labels imported. PreFilter Feature matrix shape:")
    print(X.shape)

    selectK = SelectKBest(k=kbest)
    selectK.fit(X,y)
    selectK_mask=selectK.get_support()
    K_featnames = featureNames[selectK_mask]
    print('X After K filter:',X.shape)
    print("K_featnames: %s" %(K_featnames))
    if reduceMatrix ==True :
        Reduced_df = pd.read_csv(filename, index_col=0)
        Reduced_df = Reduced_df[Reduced_df.columns[selectK_mask]]
        Reduced_df.to_csv('REDUCED_Feat.csv')
        print('Saved to REDUCED_Feat.csv')
        return Reduced_df

#WORKS! But unreadable with too many features!
def PlotFeaturesImportance(X,y,featureNames):
    '''

    Plot the relative contribution/importance of the features.
    Best to reduce to top X features first - for interpretability
    Code example from:
    http://bugra.github.io/work/notes/2014-11-22/an-introduction-to-supervised-learning-scikit-learn/
    '''

    gbc = GradientBoostingClassifier(n_estimators=100)
    gbc.fit(X, y)
    # Get Feature Importance from the classifier
    feature_importance = gbc.feature_importances_
    # Normalize The Features
    feature_importance = 100.0 * (feature_importance / feature_importance.max())
    sorted_idx = np.argsort(feature_importance)
    pos = np.arange(sorted_idx.shape[0]) + .5
    plt.figure(figsize=(16, 12))
    plt.barh(pos, feature_importance[sorted_idx], align='center', color='#7A68A6')

    #plt.yticks(pos, np.asanyarray(df.columns.tolist())[sorted_idx]) #ORIG
    plt.yticks(pos, np.asanyarray(featureNames)[sorted_idx])

    plt.xlabel('Relative Importance')
    plt.title('Features Importance')
    plt.show()

'Works!'
def PlotPerfPercentFeatures(X,y,est=LinearSVC()):
    '''
    Performance of a classifier (default: SVM-Anova)
    varying the percentile of features selected (F-test) .

    http://scikit-learn.org/stable/auto_examples/svm/plot_svm_anova.html#example-svm-plot-svm-anova-py

    See Also: (Similar but with model seelction from among classifiers):
    http://nbviewer.ipython.org/github/bugra/pydata-nyc-2014/blob/master/6.%20Scikit%20Learn%20-%20Model%20Selection.ipynb

    '''
    transform = SelectPercentile(f_classif)

    clf = Pipeline([('anova', transform), ('est', est)])
    ###############################################################################
    # Plot the cross-validation score as a function of percentile of features
    score_means = list()
    score_stds = list()
    percentiles = (1,2,3,5,7,10,13,15,20,25,33,50,65,75,90, 99)
    # percentiles = (1,5,10,25,50,75,90)

    for percentile in percentiles:
        # print(percentile)
        clf.set_params(anova__percentile=percentile)
        this_scores = cross_val_score(clf, X, y,cv=StratifiedShuffleSplit(y, n_iter=7, test_size=0.3), n_jobs=-1)
        score_means.append(this_scores.mean())
        score_stds.append(this_scores.std())
    print("Outputting Graph:")

    plt.errorbar(percentiles, score_means, np.array(score_stds))

    plt.title(
        'Predictor Performance, varying percent of features used')
    plt.xlabel('Percentile')
    plt.ylabel('Prediction Performance')
    plt.axis('tight')
    plt.show()


'NOT tested'
def plotRFECV (X,y,stepSize=0.05,scoring='f1'):
    '''
    Plot recursive feature elimination example with automatic tuning of the number of features selected with cross-validation.
    http://scikit-learn.org/stable/auto_examples/plot_rfe_with_cross_validation.html#example-plot-rfe-with-cross-validation-py
    '''
    from sklearn.svm import SVC
    from sklearn.cross_validation import StratifiedKFold
    from sklearn.feature_selection import RFECV

    # Create the RFE object and compute a cross-validated score.
    # svc = SVC(kernel="linear")
    svc = SVC(kernel="linear",class_weight='auto', cache_size=1400)
    # The "accuracy" scoring is proportional to the number of correct
    # classifications
    rfecv = RFECV(estimator=svc, step=stepSize, cv=StratifiedKFold(y, 2),
                  scoring=scoring)
    rfecv.fit(X, y)

    print("Optimal number of features : %d" % rfecv.n_features_)

    # Plot number of features VS. cross-validation scores
    import matplotlib.pyplot as plt
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Cross validation score (nb of correct classifications)")
    plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
    plt.show()

    return rfecv


def report(grid_scores, n_top=2) :
    '''
    Print out top models/parameters after a grid search for model params.
    '''
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores) :
        if n_top>1:
            print("Model with rank: {0}".format(i + 1))
        print("Average Cross-Validation score (while tuning): {0:.2f} (std: {1:.2f})".format(
            score.mean_validation_score, np.std(score.cv_validation_scores)))
        print("Model Parameters: {0}".format(score.parameters))
        print("")


def ModelParam_GridSearch(X_train, y_train, cv=4,scoreParam = 'f1'):
    '''
    Basic grid searchCV for multiple classifiers' perf & parameters.
    This is very limited and computationally expensive.
    Not guaranteed to reach even a local optima, but good to get a
    rough idea of parameters for the classifiers. (Does not address pre-processing)
    More classifiers can be added as desired, and parameters expanded.

    Later: Add options for RBM + Logit; PCA; ICA; LDA.
     See also
    http://scikit-learn-laboratory.readthedocs.org/en/latest/_modules/skll/learner.html

    TODO: Add parameters + put classifiers/"pipeline_#" in a list. (To allow checking only some params)
    '''

#    pipeline1 = Pipeline('clf', RandomForestClassifier() )
#
#    pipeline2 = Pipeline(
#    ('clf', KNeighborsClassifier()),)
    pipeline1 = RandomForestClassifier(n_jobs=-1)
    pipeline2 = KNeighborsClassifier()

    pipeline3 = SVC(cache_size=1500)
    # pipeline3 = NuSVC(cache_size=1500)

    pipeline4 = GaussianNB()
    pipeline5 = GradientBoostingClassifier()
    pipeline6 = SGDClassifier()
    pipeline7 = LogisticRegression()


    'RandomForestClassifier:'
    parameters1 = {
    'n_estimators': [150],
    'criterion': ['gini'],
    'max_features': ['auto',0.4],
    'max_depth': [8,None],
    'min_samples_leaf':[1,2],
    'min_samples_split':[2,4],
    'n_jobs':[-1]
    }
    #, 'entropy'
        # 'n_jobs':[-1]

    'KNeighborsClassifier:'
    parameters2 = {
    'n_neighbors': [7],
    'weights': ['distance']
    }

    'SVC:'
    parameters3 = {
    'C': [0.01,0.1, 1,10,100],
    'kernel': ['linear','rbf'],
    'gamma': [0.1,0.0, 1.0,20],
    'cache_size':[1500],
    'class_weight':['auto'],
    }
# , 'poly','sigmoid']

##    'GaussianNB:'
##    parameters4 = {}

    'GradientBoostingClassifier'
    parameters5 = {
    'max_depth':[3,5,8],
    'n_estimators': [100],
    'min_samples_leaf':[1,2],
    'learning_rate': [0.1, 0.01],
    'max_features': ['auto',0.4]
    }
    'SGDClassifier:'
    parameters6 = {
     'alpha': [0.00001,0.001,0.01],
    'penalty': ['l1','l2', 'elasticnet'],
    'n_iter': [300],
    'loss':['hinge'],
    'n_jobs':[-1],
    'class_weight':['auto']
    }
#, 'modified_huber','log'

    'LogisticRegression:'
    parameters7 = {
    'C': [0.001,0.01, 0.1, 1.0,10,100],
    'penalty': ['l1','l2'],
    'class_weight':['auto']
    }

    'TODO: make this into a seperate method, with pars, pips passed to it as params'
    pars = [parameters1, parameters2, parameters3,parameters5,parameters6,parameters7] #parameters4
    pips = [pipeline1, pipeline2, pipeline3,pipeline5,pipeline6,pipeline7] # pipeline4,

    print ("Starting Gridsearch To find each model's best parameters")
    for i in range(len(pars)):
        print(pips[i])

        gs = GridSearchCV(estimator=pips[i], param_grid=pars[i],
                          verbose=0, refit=True, n_jobs=-1,iid=False,
                          pre_dispatch='2*n_jobs',scoring=scoreParam,
                          fit_params={'sample_weight': balance_weights(y)},
                          cv=StratifiedKFold(y_train,n_folds=cv,shuffle=True))
#Valid scoring options: ['accuracy', 'average_precision', 'f1', 'precision', 'recall', 'roc_auc']
        # gs = gs.fit(X_train, y_train)
        'http://stackoverflow.com/questions/13051706/scikit-learn-using-sample-weight-in-grid-search?rq=1'
        'Note: Remove "class_weight=auto"  from the autoweighting classifiers!!'
        "Set Class weights (then into sample weights: https://github.com/scikit-learn/scikit-learn/blob/8dab222cfe894126dfb67832da2f4e871b87bce7/sklearn/utils/class_weight.py"
        gs.fit(X_train, y_train)
        #print ("Finished Gridsearch")
        #print (gs.best_score_)
        report(gs.grid_scores_)
        # http://stackoverflow.com/questions/18210799/scikit-learn-sample-try-out-with-my-classifier-and-data

        'Get more exhaustive CV results with the best tuned parameters for the model'
        est = gs.best_estimator_
        scores = cross_val_score(est, X_train,
                                 y_train,
                                 cv=StratifiedShuffleSplit(y=y_train, n_iter=10, test_size=0.2),scoring=scoreParam,
                                 n_jobs=-1, pre_dispatch='1.8*n_jobs')
        print("Tuned Model's %s Score: %0.3f (+/- %0.3f)" % (scoreParam,scores.mean(), scores.std() * 2))



if __name__ == '__main__' :
    'TODO: Allow user to select desired function - CV model, or feature reduction'
    'TODO: Use os.path.join - for file names/locations/dirs..'
    #Set by user input:
    fileName = r'/trainingSetFeatures.csv'
    filePath = str(argv[1])
    X, y, lb_encoder,featureNames = load_data(filePath+fileName, 'file') # X, y = features, labels

    print(X.shape,"= (samples, features)")

    y_inv = Counter(lb_encoder.inverse_transform(y))
    print("Classes:", y_inv)

    # 'Normalize/Scale features if needed. Our data is standardized by default'
    # X = StandardScaler(copy=False).fit_transform(X)

    Fwe = SelectFwe(alpha=0.01).fit(X,y)
    X=Fwe.transform(X)
    featureNames=featureNames[Fwe.get_support()]
    print("F-test filter ->",X.shape)

    FeatSelection_SVM=True
    FeatSelection_RandLogReg=False

    if FeatSelection_RandLogReg == True:
        LogRegFeats = RandomizedLogisticRegression(C=5, scaling=0.5,
         sample_fraction=0.8, n_resampling=60, selection_threshold=0.2,n_jobs=-1)
        X = LogRegFeats.fit_transform(X,y)
        featureNames=featureNames[LogRegFeats.get_support()]
        print("RandomizedLogisticRegression Feature Selection ->:",X.shape)

    elif FeatSelection_SVM == True:
        X= LinearSVC(C=1, penalty="l1", dual=False,class_weight='auto').fit_transform(X, y)
        # X= LogisticRegression(C=0.01,class_weight='auto').fit_transform(X, y)
        featureNames=featureNames[LogRegFeats.get_support()]
        print ("SVC Transformed X:",X.shape)

    '''
    print("Plot #Feats vs Classification performance:")
    PlotPerfPercentFeatures(X_LR,y,est=SVC(C=100))
    '''

    KFilt=None
    # KFilt=200

    if KFilt is not None:
        k = SelectKBest(k=KFilt).fit(X,y)
        X=k.transform(X)
        featureNames=featureNames[k.get_support()]
        print("X reduced to K best features: ",X.shape)

    print("Performance as a function of percent of features used:")
    PlotPerfPercentFeatures(X,y,est=LinearSVC())

    #varFilt = VarianceThreshold(threshold=0.05)
    #X = varFilt.fit_transform(X)
    #print(X.shape,"X post low variance feature filtering")

    'EG - graph best features; feature selection using RF, ensemble classifiers..'
    'http://nbviewer.ipython.org/github/herrfz/dataanalysis/blob/master/assignment2/samsung_data_prediction_submitted.ipynb'

    RFE_FeatsToKeep = 15
    FeatSelection_RFE=True
    FeatSelection_RFECV=False

    if (FeatSelection_RFE or FeatSelection_RFECV) == True:
        'RFE + - best feats'
        'http://scikit-learn.org/stable/auto_examples/plot_rfe_with_cross_validation.html '
        svc = LinearSVC(class_weight='auto')#,penalty='l1',dual=False)
        # svc = LogisticRegression(class_weight='auto')#,C=1)
        if FeatSelection_RFECV==True:
            rfecv = RFECV(estimator=svc, step=0.1,
                         cv=StratifiedShuffleSplit(y,n_iter=7,test_size=0.33),
                         scoring='f1',verbose=0)
            # " scoring='roc_auc','recall','f1'..."
        else:
            rfecv = RFE(estimator=svc,n_features_to_select=RFE_FeatsToKeep, step=0.1)
        rfecv.fit(X, y)
        if FeatSelection_RFECV==True:
            print("RFEcv selected %d number of Optimal features : " % (rfecv.n_features_))
        print("RFE (%d Features) scorer : \n" % (rfecv.n_features_),rfecv.score(X, y) )
        print("RFE selected feature names:")
        featureNames=featureNames[rfecv.get_support()]
        rfe_featnames = featureNames[rfecv.get_support()]
        print (rfe_featnames)
        X_RFE = rfecv.fit_transform(X, y)
        print(X_RFE.shape,"X_RFE \n")

        'Set GetRFEPerf To true or by user, if perf. of reduced set wanted'
        GetRFEPerf=False


    print("\n X: \n")
    ModelParam_GridSearch(X,y,cv=4)

    if GetRFEPerf==True:
        print("\n X-RFE: \n")
        ModelParam_GridSearch(X_RFE,y,cv=4)

    GetPCAPerf=False
    if GetPCAPerf==True:
        pca = PCA(n_components=0.99,whiten=False)
        X_PCA = pca.fit_transform(X)
        print(X_PCA.shape,"X + PCA")
        print("X_PCA \n")
        ModelParam_GridSearch(X_PCA,y,cv=3)
