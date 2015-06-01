#!/cs/prt3/danofer/anaconda3/bin/python3
# coding: utf-8

# "http://nbviewer.ipython.org/github/gmonce/scikit-learn-book/blob/master/Chapter%204%20-%20Advanced%20Features%20-%20Feature%20Engineering%20and%20Selection.ipynb"

import os
import sys
import argparse
import pandas
import numpy
import operator

# import seaborn only on systems which have it installed
try:
    import seaborn as sns
except ImportError:
    pass

import matplotlib as mpl
import matplotlib.pyplot as plt

from sklearn.preprocessing import LabelEncoder
from sklearn import cross_validation
from sklearn.feature_selection import RFE, RFECV,SelectFwe, SelectKBest
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.svm import LinearSVC, SVC

class RandomForestClassifierWithCoef(RandomForestClassifier):
    '''
    Let's us use RF with RFE feature selection.
    http://stackoverflow.com/questions/24123498/recursive-feature-elimination-on-random-forest-using-scikit-learn?rq=1
    '''

    def fit(self, *args, **kwargs):
        super(RandomForestClassifierWithCoef, self).fit(*args, **kwargs)
        self.coef_ = self.feature_importances_


def report(grid_scores, n_top=1) :
    '''
    Print out top models/parameters after a grid search for model params.
    '''
    top_scores = sorted(grid_scores, key=operator.itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores) :
        if n_top>1:
            print("Model with rank: {0}".format(i + 1))
        print("Average Cross-Validation score (while tuning): {0:.3f} (std: {1:.2f})".format(
            score.mean_validation_score, numpy.std(score.cv_validation_scores)))
        print("Model Parameters: {0}".format(score.parameters))
        print("")

def PlotFeaturesImportance(X,y,featureNames,dataName):
    '''
    Plot the relative contribution/importance of the features.
    Best to reduce to top X features first - for interpretability
    Code example from:
    http://bugra.github.io/work/notes/2014-11-22/an-introduction-to-supervised-learning-scikit-learn/
    '''
    gbc = GradientBoostingClassifier(n_estimators=40)
    gbc.fit(X, y)
    # Get Feature Importance from the classifier
    feature_importance = gbc.feature_importances_
    # Normalize The Features
    feature_importance = 100 * (feature_importance / feature_importance.max())
    sorted_idx = numpy.argsort(feature_importance)
    pos = numpy.arange(sorted_idx.shape[0]) + 4.5
    # pos = numpy.arange(sorted_idx.shape[0])
    # plt.figure(figsize=(16, 12))
    plt.figure(figsize=(14, 9), dpi=250)
    plt.barh(pos, feature_importance[sorted_idx], align='center', color='#7A68A6')
    #plt.yticks(pos, numpy.asanyarray(df.columns.tolist())[sorted_idx]) #ORIG
    plt.yticks(pos, numpy.asanyarray(featureNames)[sorted_idx])

    plt.xlabel('Relative Importance')
    plt.title('%s: Top Features' %(dataName))
    plt.grid('off')
    plt.ion()
    plt.show()
    plt.savefig(str(dataName)+'TopFeatures.png',dpi=200)

def altPlotFeaturesImportance(X,y,featureNames,dataName):
    "http://nbviewer.ipython.org/github/cs109/2014/blob/master/homework-solutions/HW5-solutions.ipynb"
    clf = RandomForestClassifier(n_estimators=50)

    clf.fit(X,y)
    importance_list = clf.feature_importances_
    # name_list = df.columns #ORIG
    name_list=featureNames

    importance_list, name_list = zip(*sorted(zip(importance_list, name_list)))
    plt.barh(range(len(name_list)),importance_list,align='center')
    plt.yticks(range(len(name_list)),name_list)
    plt.xlabel('Relative Importance in the Random Forest')
    plt.ylabel('Features')
    plt.title('%s \n Relative Feature Importance' %(dataName))
    plt.grid('off')
    plt.ion()
    plt.show()

def PlotPerfPercentFeatures(X,y,est=LinearSVC()):
    '''
    Display performance of a classifier (default: SVM),
    varying the percentile of features retained (F-test) .

    http://scikit-learn.org/stable/auto_examples/svm/plot_svm_anova.html#example-svm-plot-svm-anova-py
    '''
    transform = SelectPercentile(f_classif)

    clf = Pipeline([('anova', transform), ('est', est)])
    ###############################################################################
    # Plot the cross-validation score as a function of percentile of features
    score_means = list()
    score_stds = list()
    percentiles = (1,2,3,5,7,10,15,20,25,33,50,66,75,90, 100)
    # percentiles = (1,5,10,25,50,75,90)

    for percentile in percentiles:
        # print(percentile)
        clf.set_params(anova__percentile=percentile)
        this_scores = cross_val_score(clf, X, y,cv=StratifiedShuffleSplit(y, n_iter=5, test_size=0.4), n_jobs=-1)
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

'''
# http://stackoverflow.com/questions/27264542/get-row-and-column-names-argmax-for-max-entry-in-pandas-dataframe
# Remove correlated features

S = S.abs()
np.fill_diagonal(S.values,-1) # so that max can't be on the diagonal now
S = rec_drop(S,max_allowed_correlation=0.95)

def rec_drop(S, max_allowed_correlation=0.99):
    max_corr = S.max().max()
    if max_corr<max_allowed_correlation: # base case for recursion
         return S.columns.tolist()
    row,col = arg_max(S)  # row and col are distinct features - max can't be on the diagonal
    S = S.drop(row).drop(row,axis=1) # removing one of the features from S
    return rec_drop(S, max_allowed_correlation)
'''

def PlotConfusionMatrix (cm,target_names,Title="Neuropeptides"):
    '''
    Works!
    '''
    import matplotlib as mpl
    np.set_printoptions(suppress=True)
    mpl.rc("figure", figsize=(12, 9))

    hm = sns.heatmap(cm,
                cbar=True,
                annot=True,
                square=True,
                fmt='d',
                yticklabels=target_names,
                xticklabels=target_names,
                # cmap='Blues'
                )
    plt.fig(dpi=200)
    plt.title('Confusion matrix: '+Title)
    plt.ylabel('Actual class')
    plt.xlabel('Predicted class')
    plt.tight_layout()
##    plt.savefig('./images/confmat.png', dpi=300)
    plt.show()


def main(args):
    if args.train_dir is None:
        # args.train_dir = '/a/fr-05/vol/protein/danofer/ProtFeat/feat_extract/chap/train/'
        #args.train_dir = '/cs/prt3/danofer/ProtFeat/feat_extract/test_seq/NP/SPCleaved_NP-70+NEG-30_Big-V3/'
#        args.train_dir =  r'D:\SkyDrive\Dropbox\bioInf_lab\AA_info\CODE\feat_extract\test_seq\NP\SPCleaved_NP-70+NEG-30_Big-V3'
        # args.train_dir =  r'E:\Dropbox\Dropbox\bioInf_lab\AA_info\fastas\NP\SP_Cleaved+NP+Neg_Big'
        args.train_dir =  r'E:\Dropbox\Dropbox\bioInf_lab\AA_info\fastas\Benchmarks\Thermophiles'
        print("Using default train_dir: %s" % args.train_dir)

    pandas.set_option('display.max_columns', 10)
    pandas.set_option('display.max_rows', 4)
    # mpl.rc('title', labelsize=6)
    mpl.rc('ytick', labelsize=7)
    mpl.rc('xtick', labelsize=4)

    os.chdir(args.train_dir)
    dataName = 'Neuropeptides'

    df = pandas.read_csv('trainingSetFeatures.csv')
    feature_cols = [col for col in df.columns if col not in ['classname','Id','proteinname']]
    feature_cols=numpy.array(feature_cols)

    X = df[feature_cols].values
    y = df.classname.values

    le = LabelEncoder()
    y = le.fit_transform(y)

    "Initial feature selection trimming"
    print(X.shape)

    Fwe = SelectFwe(alpha=0.01).fit(X,y)
    X=Fwe.transform(X)
    print("F-test -> ",X.shape)
    feature_cols=feature_cols[Fwe.get_support()]
    '''
    FeatSelection_SVM = True
    if FeatSelection_SVM == True:
        svc_L1 = LinearSVC(C=50, penalty="l1", dual=False,class_weight='auto').fit(X, y)
        X = svc_L1.transform(X, y)
        print ("L1 SVM Transformed X:",X_L1.shape)
        feature_cols=feature_cols[list(set(np.where(svc_L1.coef_ != 0)[-1]))]
    '''


    k = SelectKBest(k=255).fit(X,y)
    X=k.transform(X)
    feature_cols=feature_cols[k.get_support()]


    param_dist = {"max_depth": [6,9, None],
                  "max_features": ['auto',0.4],
                  "min_samples_leaf": [1,2,3],
                  "bootstrap": [True, False],
                  'min_samples_split':[2,3],
                  "criterion": [ "gini"],
                  "n_estimators":[100],
                  "n_jobs":[-1]}

    rf = RandomForestClassifierWithCoef(max_depth= 7, min_samples_split= 1, min_samples_leaf= 2, n_estimators= 50,  n_jobs= 2, max_features= "auto")

    "WARNING! F1 Score as implemented by Default in binary classification (two classes) gives the score for 1 class."

    scores = cross_validation.cross_val_score(rf,X,y,n_jobs=-1,cv=cross_validation.StratifiedShuffleSplit(y,n_iter=8,test_size=0.2))
    print("X RF Accuracy: %0.3f (+- %0.2f)" % (scores.mean(), scores.std() * 2))
    "Instead of scores_f1, we could also use precision, sensitivity, MCC (if binary), etc'."
    scores_f1 = cross_validation.cross_val_score(rf,X,y,n_jobs=-1,cv=cross_validation.StratifiedShuffleSplit(y,n_iter=8,test_size=0.2),scoring='f1')
    print("X RF f1: %0.3f (+- %0.2f)" % (scores_f1.mean(), scores_f1.std() * 2))

    # rfeSelect = RFE(estimator=rf,n_features_to_select=16, step=0.04)
    rfeSelect = RFECV(estimator=rf,step=20, cv=2,scoring='f1') #average_precision , recall
    X_RFE = rfeSelect.fit_transform(X,y)
    print(X_RFE.shape)

    RFE_FeatureNames = feature_cols[rfeSelect.get_support()]
    print(RFE_FeatureNames)

    RFE_ScoreRatio = 100*(cross_validation.cross_val_score(rf,X_RFE,y,n_jobs=-1,cv=cross_validation.StratifiedShuffleSplit(y,n_iter=8,test_size=0.2),scoring='f1').mean())/scores_f1.mean()
    print("Even with just",X_RFE.shape[1]," features, we have %f performance! (f1 score ratio)" %(RFE_ScoreRatio))

    # PlotFeaturesImportance(X_RFE, y, RFE_FeatureNames, dataName)
    print("Alt plot:")
    altPlotFeaturesImportance(X_RFE, y, RFE_FeatureNames, dataName)

def read_cmd(args):
    ap = argparse.ArgumentParser()
    ap.add_argument('-t', '--train_dir', help="Set train dir",
                    default=None)
    return ap.parse_args(args)


if __name__ == '__main__':
    sys.exit(main(read_cmd(sys.argv[1:])))
