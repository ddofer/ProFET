# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 18:24:07 2014

@author: Dan
Modified cut from PipeTasks - tries only 2-3 models and outputs results + dummy res. + stats. (+ foldername)
Common tasks in processing data and models, likely to be called as part of the pipeline.
"""

#Lab workaround
import sys
sys.path += [r'C:\Anaconda\Lib', r'C:\Anaconda\Lib\site-packages']
sys.path += [r'/cs/prt3/danofer/anaconda3',r'/cs/prt3/danofer/anaconda3/Lib']


from sklearn import metrics
from sys import argv
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier, ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression,  RandomizedLogisticRegression
from sklearn.svm import SVC,LinearSVC
from sklearn.naive_bayes import GaussianNB
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import StratifiedKFold,cross_val_score,StratifiedShuffleSplit
# from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE, RFECV, SelectFdr,SelectFwe,SelectKBest
from sklearn.linear_model import RandomizedLogisticRegression
from sklearn.lda import LDA
from sklearn.decomposition import PCA#,FastICA #,TruncatedSVD
from operator import itemgetter
from Model_trainer import load_data
from sklearn.metrics import precision_recall_fscore_support,matthews_corrcoef, classification_report,confusion_matrix
#http://tokestermw.github.io/posts/imbalanced-datasets-random-forests/
# from sklearn.preprocessing import balance_weights
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.metrics import f1_score
from PipeTasks import Get_yPred,balance_weights
from collections import Counter
from sklearn.dummy import DummyClassifier
import os, fnmatch
import seaborn as sns
from sklearn.preprocessing import LabelEncoder

'Use SKLL!? https://skll.readthedocs.org/en/latest/run_experiment.html#quick-example'
'''
 enhanced confusion matrix + weights for RF:
http://stackoverflow.com/questions/24123498/recursive-feature-elimination-on-random-forest-using-scikit-learn
'''

'Works!'
def PlotPerfPercentFeatures(X,y,est=LinearSVC()):
    '''
    Performance of a classifier (default: SVM-Anova)
    varying the percentile of features selected (F-test) .

    http://scikit-learn.org/stable/auto_examples/svm/plot_svm_anova.html#example-svm-plot-svm-anova-py
    '''
    transform = SelectPercentile(f_classif)

    clf = Pipeline([('anova', transform), ('est', est)])
    ###############################################################################
    # Plot the cross-validation score as a function of percentile of features
    score_means = list()
    score_stds = list()
    percentiles = (1,2,3,5,7,10,13,15,20,25,33,50,65,75,90, 100)
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


def PlotConfusionMatrix (cm,target_names,Title="Neuropeptides"):
    '''
    Works!
    '''
    import matplotlib as mpl
    np.set_printoptions(suppress=True)
    mpl.rc("figure", figsize=(16, 12))

    hm = sns.heatmap(cm,
                cbar=True,
                annot=True,
                square=True,
                fmt='d',
                yticklabels=target_names,
                xticklabels=target_names,
                cmap='Blues'
                )
    plt.title('Confusion matrix'+Title)
    plt.ylabel('Actual class')
    plt.xlabel('Predicted class')
    plt.tight_layout()
##    plt.savefig('./images/confmat.png', dpi=300)
    plt.show()


def CV_multi_stats(X, y, model,n=6) :
    '''
    http://scikit-learn.org/stable/modules/model_evaluation.html#classification-metrics
    This version uses multiclass (or multilabel) compatible metrics.

    May be expanded to use the cross_val_score helper function:
    http://scikit-learn.org/stable/modules/generated/sklearn.cross_validation.cross_val_score.html
    http://scikit-learn.org/stable/modules/cross_validation.html#computing-cross-validated-metrics
    '''

    scores = cross_val_score(estimator=model, X=X, y=y, cv=StratifiedShuffleSplit(y, n_iter=n, test_size=0.16), n_jobs=-1) #Accuracy
    scores_f1 = cross_val_score(estimator=model, X=X, y=y, cv=StratifiedShuffleSplit(y, n_iter=n, test_size=0.16), n_jobs=-1, scoring='f1')
    print("Model Accuracy: %0.3f (+- %0.2f)" % (scores.mean(), scores.std() * 2))
    print("Model f1: %0.3f (+- %0.2f)" % (scores_f1.mean(), scores_f1.std() * 2))
    return (scores.mean(), scores.std() ,scores_f1.mean(), scores_f1.std() ) #Removed * 2 from returned STD .. ?


def CV_Binary_stats(X, y, model,n=10) :
    '''
    http://scikit-learn.org/stable/modules/model_evaluation.html#classification-metrics
    Note that some of the metrics here ONLY work for BINARY tasks.
    This will be VERY slow compared to the built-in, multicore CV implementation. (Unless
     used with a classifier that is parallelized anyway, such as RF).
    By default, balances weights when fitting

    http://scikit-learn.org/stable/modules/cross_validation.html#computing-cross-validated-metrics
    '''
    from sklearn.metrics import precision_score, accuracy_score, recall_score,precision_recall_fscore_support

    mean_auc = 0.0
    mean_precision = 0.0
    mean_recall = 0.0
    mean_accuracy = 0.0

    sss = StratifiedShuffleSplit(y,  n_iter=n, test_size=0.2, random_state=0)
    for train_index, test_index in sss:
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

    # for i in range(n) :
    #     # for each iteration, randomly hold out 30% of the data as CV set
    #     X_train, X_cv, y_train, y_cv = cross_validation.train_test_split(X, y,
    #                                                                      test_size=.15,
    #                                                                      random_state=i)
    #     cv=StratifiedShuffleSplit(y=y_train, n_iter=11, test_size=0.11)
        # train model and make predictions
        model.fit(X_train, y_train,sample_weight=balance_weights(y_train))
        # preds = model.predict(X_cv)
        preds = model.predict(X_test)

        '''
        # ROC_AUC - Restricted to binary (not multiclass) case.
        fpr, tpr, thresholds = metrics.roc_curve(y_cv, preds)
        roc_auc = metrics.auc(fpr, tpr)
        # print("( %d/%d)" % (i + 1, n))
        mean_auc += roc_auc
        '''
        accuracy = accuracy_score(y_cv, preds)
        precision = precision_score(y_cv, preds)
        recall = recall_score(y_cv, preds)
        mean_accuracy += accuracy
        mean_precision += precision
        mean_recall += recall

    mean_accuracy = (mean_accuracy / n)
    mean_precision = mean_precision / n
    mean_recall = mean_recall / n
    # mean_auc = mean_auc / n
    print('mean_accuracy:  %s ' %(round(mean_accuracy, 3)))
    print('mean_precision:  %s ' %(round(mean_precision, 3)))
    print('mean_recall:  %s ' %(round(mean_recall, 3)))
    # print('mean_auc:  %s ' %(round(mean_auc, 3)))
    return (mean_accuracy,mean_precision,mean_recall)


def find_files(directory, pattern):
    'http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python'
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


def report(grid_scores, n_top=1) :
    '''
    Print out top models/parameters after a grid search for model params.
    '''
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores) :
        if n_top>1:
            print("Model with rank: {0}".format(i + 1))
        print("Average (tuning) CV score : {0:.2f} (std: {1:.2f})".format(
            score.mean_validation_score, np.std(score.cv_validation_scores)))
        print("Optimal Model Hyperparameters: {0}".format(score.parameters))
        print("")


def ModelParam_GridSearch(X_train, y_train, cv=3,scoreParam = 'precision'):
    '''
    Basic grid searchCV for multiple classifiers' perf & parameters.
    This is limited as currently implemented, but still  computationally expensive.
    Not guaranteed to reach even a local optima, but good to get a
    rough idea of parameters for the classifiers. (Does not address pre-processing)
    More classifiers can be added as desired, and parameters expanded.

    Later: Add options for RBM + Logit; PCA; ICA; LDA.
    (Further) Feature selection should  be implemented within the CV pipeline, if you wish
    to avoid overfitting. (Note the effects of )
     See also
    http://scikit-learn-laboratory.readthedocs.org/en/latest/_modules/skll/learner.html

    Possible Scoreparams: scoreParam = 'f1','accuracy', 'precision', 'roc_auc'..
    '''

#    pipeline1 = Pipeline('clf', RandomForestClassifier() )

    pipeline1 = RandomForestClassifier(n_jobs=-1)
    pipeline2 = SVC(cache_size=1900)
    pipeline3 = GradientBoostingClassifier()
    pipeline4 = LogisticRegression()


    'RandomForestClassifier:'
    parameters1 = {
    'n_estimators': [120],
    'criterion': ['gini'],
    'max_features': ['auto',0.4],
    'min_samples_leaf':[1,2],
   'min_samples_split':[2,3],
    'n_jobs':[-1],
    'max_depth': [8, None]

    }

    'SVC:'
    parameters2 = {
    'C': [0.2, 1,10,50,100,1000],
    # 'kernel': ['linear','rbf'],
    'kernel': ['rbf'],
    'gamma': [0.1,0.0, 1.0],
    'cache_size':[1900],
    'class_weight':['auto',None],
    }
# , 'poly','sigmoid']
    'GradientBoostingClassifier'
    parameters3 = {
    'max_depth':[5,7],
    'n_estimators': [80],
    # 'min_samples_leaf':[2],
    # 'learning_rate': [0.1, 0.05],
    'max_features': ['auto',0.4]
    }
    # 'min_samples_leaf':[1,2],

    'LogisticRegression:'
    parameters4 = {
    'C': [1.0,10,100],
    'penalty': ['l1','l2'],'class_weight':['auto',None]
    }

    pars = [parameters1, parameters2, parameters3,parameters4]
    pips = [pipeline1, pipeline2, pipeline3, pipeline4]

    'Store and return the best estimator found (and score)'
    bestEst=None
    bestScore=0

    print ("Starting gridsearch to find best model hyperparameters.")

    'Gridsearch done "in bits" due to some classifiers not supporting sample_Weight'
    def gs_fit(gs):
        nonlocal bestEst
        nonlocal bestScore
        gs.fit(X_train, y_train)
        report(gs.grid_scores_)
        # http://stackoverflow.com/questions/18210799/scikit-learn-sample-try-out-with-my-classifier-and-data
        if gs.best_score_>bestScore:
            bestEst = gs.best_estimator_
            bestScore = gs.best_score_
            print("Updated best Est, new Best score:",bestScore)

    # for i in range(len(pars)): #Orig
    for i in range(2):
        clf_name = str(pips[i])
        print(clf_name[0:clf_name.index("(")])
        gs = GridSearchCV(estimator=pips[i], param_grid=pars[i],
                          verbose=1, refit=True, n_jobs=-1,iid=False,
                          fit_params={'sample_weight': balance_weights(y_train)},
                          pre_dispatch='1.5*n_jobs',scoring=scoreParam,
                          cv=StratifiedKFold(y_train,n_folds=cv,shuffle=True))
#Valid scoring options: ['accuracy', 'average_precision', 'f1', 'precision', 'recall', 'roc_auc']
        gs_fit(gs)

    i=3 #Logistic Regression

    gs = GridSearchCV(estimator=pips[i], param_grid=pars[i],
                          verbose=0, refit=True, n_jobs=-1,iid=True,
                          pre_dispatch='1.5*n_jobs',scoring=scoreParam,
                          cv=StratifiedKFold(y_train,n_folds=cv,shuffle=True))
    gs_fit(gs)

   # 'http://stackoverflow.com/questions/13051706/scikit-learn-using-sample-weight-in-grid-search?rq=1'
    # "http://stackoverflow.com/questions/20082674/unbalanced-classification-using-randomforestclassifier-in-sklearn"
    # "Set Class weights (then into sample weights: https://github.com/scikit-learn/scikit-learn/blob/8dab222cfe894126dfb67832da2f4e871b87bce7/sklearn/utils/class_weight.py"
    #print (gs.best_score_)

    print("Best Predictor:", bestEst, "Score: ",(bestScore))
    return(bestEst,bestScore)

def fileNameFromPaths (filePaths):
    '''
    Warning: VERY Ad-hoc way of getting the name without the filepath. Will breake easily
    '''
    names =[]
    for filePath in filePaths:
        filePath = os.path.normpath(filePath)
        fileName = str(filePath)
        # fileName = fileName[fileName.rindex(r"test_seq/"):fileName.index('trainingSetFeatures')]
        names.append(fileName[fileName.find("test_seq/")+9:fileName.rindex('trainingSetFeatures')-1])
    return names

def GetAllPerf (filePaths=None):
    if filePaths is None:
        filePaths = list(find_files(directory='./test_seq', pattern='trainingSetFeatures.csv'))

    #Sanity check:
    # filePaths=['/a/fr-05/vol/protein/danofer/ProtFeat/feat_extract/test_seq/Thermophile']
    # filePaths=['./test_seq/NP/NP2/Train/trainingSetFeatures.csv']

    print("FilePaths: \n",filePaths)
    fileNames=fileNameFromPaths (filePaths)
    print("FileNames:",fileNames)


    resDict = pd.DataFrame(index=fileNames,
        columns=['Accuracy','Accuracy_SD',
        'f1','f1_SD','dummy_freq:Accuracy','dummy_freq:f1',
        'LargestClassPercent','Classes',
        # 'TopRFE-Features','Best (f1) Model parameters',
         '# Classes',
         'Array-Acc-Scores' ,'Array-f1-Scores'
         ,'bestML-Acc','bestML-f1','dummy_freq_f1_weighted'])


    #redDict holds results for each file/class, for saving to output-file

    i=-1
    for filePath in filePaths:
        i +=1

        'http://pythonconquerstheuniverse.wordpress.com/2008/06/04/gotcha-%E2%80%94-backslashes-in-windows-filenames/'
        filePath = os.path.normpath(filePath)
        print(filePath)
        fileName=str(fileNames[i]) #Str added now 14.1

        print("fileName: %s" %(fileName))
        "resDict['Name']= fileName"

        # filePath = str(argv[1])
        # X, y, lb_encoder,featureNames = load_data(filePath+fileName, 'file') # X, y = features, labels
        X, y, lb_encoder,featureNames = load_data(filePath) # X, y = features, labels
        print(X.shape,"= (samples, features)")
        y_inv = Counter(lb_encoder.inverse_transform(y))
        MajorityPercent = round(100*y_inv.most_common()[0][1]/sum(y_inv.values()),1)
        print("Classes:", lb_encoder.classes_)
        print("MajorityClassPercent:", MajorityPercent)

        resDict.LargestClassPercent[fileName] = MajorityPercent
        resDict.Classes[fileName] = str(lb_encoder.classes_)
        resDict["# Classes"][fileName]=len(lb_encoder.classes_)

        KFilt=None
        KFilt=350  #This is just temporary for the outputs - saves computation time. Barely filters compared to the model itself.

        if KFilt is not None:
            k = SelectKBest(k=KFilt).fit(X,y)
            X=k.transform(X)
            featureNames=featureNames[k.get_support()]

        Fwe = SelectFwe(alpha=0.01).fit(X,y)
        X=Fwe.transform(X)
        featureNames=featureNames[Fwe.get_support()]

        print("X reduced to K best features: ",X.shape)


        FeatSelection_SVM=False #Feature Names need updating!!
        FeatSelection_RandLogReg=False

        if FeatSelection_RandLogReg == True:
            LogRegFeats = RandomizedLogisticRegression(C=10, scaling=0.5,
             sample_fraction=0.95, n_resampling=40, selection_threshold=0.2,n_jobs=-1).fit(X,y)
            X_L1 = LogRegFeats.transform(X)
            featureNames=featureNames[LogRegFeats.get_support()]
            print("RandomizedLogisticRegression Feature Selection ->:",X_L1.shape)

        elif FeatSelection_SVM == True:
            svc_L1= LinearSVC(C=30, penalty="l2", dual=False,class_weight='auto').fit(X, y)
            X_L1 = svc_L1.transform(X, y)
            featureNames=featureNames[list(set(np.where(svc_L1.coef_ != 0)[-1]))]
            print ("L1 SVM Transformed X:",X_L1.shape)
        # X=X_L1

        '''
        print("Performance as a function of percent of features used:")
        PlotPerfPercentFeatures(X,y,est=LinearSVC())
        '''

        'EG - graph best features; feature selection using RF, ensemble classifiers..'
        'http://nbviewer.ipython.org/github/herrfz/dataanalysis/blob/master/assignment2/samsung_data_prediction_submitted.ipynb'

        RFE_FeatsToKeep = 16
        FeatSelection_RFE=False
        FeatSelection_RFECV=False

        if (FeatSelection_RFE or FeatSelection_RFECV) == True:
            'RFE + - best feats'
            'http://scikit-learn.org/stable/auto_examples/plot_rfe_with_cross_validation.html '
            svc = LinearSVC(class_weight='auto')#,penalty='l1',dual=False)
            # svc = LogisticRegression(class_weight='auto')#,C=1)

            if FeatSelection_RFECV==True:
                rfecv = RFECV(estimator=svc, step=RFE_FeatsToKeep,scoring='average_precision')
                             # ,cv=StratifiedShuffleSplit(y,n_iter=3,test_size=0.3))
                             #,scoring='f1',verbose=0) # " scoring='roc_auc','recall','f1',accuracy..."
            else:
                rfecv = RFE(estimator=svc,n_features_to_select=RFE_FeatsToKeep, step=0.03)
            rfecv.fit(X, y)
            if FeatSelection_RFECV==True:
                print("RFE-CV selected %d features : " % (rfecv.n_features_))
            print("RFE (%d features) scorer : " % (rfecv.n_features_),rfecv.score(X, y) )
            rfe_featnames = featureNames[rfecv.get_support()]
            featureNames = featureNames[rfecv.get_support()]
            print("RFE selected feature names:",rfe_featnames)
            X_RFE = rfecv.fit_transform(X, y)
            print("X_RFE",X_RFE.shape)

            resDict['TopRFE-Features'][fileName]=str(rfe_featnames)

            'Set GetRFEPerf To true or by user, if perf. of reduced set wanted'
        GetRFEPerf=False

        # print("lb_encoder.classes_",lb_encoder.classes_)
        'Blind score boxplot graphic example using Seaborn: http://nbviewer.ipython.org/github/cs109/2014/blob/master/homework-solutions/HW5-solutions.ipynb '
        'Confusion matrixes + Dummies - http://bugra.github.io/work/notes/2014-11-22/an-introduction-to-supervised-learning-scikit-learn/'
        'http://scikit-learn.org/stable/modules/model_evaluation.html#dummy-estimators'

        "http://blog.yhathq.com/posts/predicting-customer-churn-with-sklearn.html"
        print()

        "Make custom F1 scorer. May not have fixed problem!"
        from sklearn.metrics.score import make_scorer
        f1_scorer = make_scorer(metrics.f1_score,
                     greater_is_better=True, average="micro") #Maybe another metric? May NOT be fixed!?. #weighted, micro, macro, none

        # print("Dummy classifiers output:")

        dummy_frequent = DummyClassifier(strategy='most_frequent',random_state=0)
        y_dummyPred = Get_yPred(X,y,clf_class=dummy_frequent)
        dummy_freq_acc = '{:.3}'.format(metrics.accuracy_score(y,y_dummyPred ))
        dummy_freq_f1 = '{:.3}'.format(metrics.f1_score(y, y_dummyPred,average='weighted'))

        dummy_freq_f1_weighted = '{:.3}'.format(f1_scorer(y, y_dummyPred))
        #Get from ALL classes f1..
        dummy_freq_f1_mean=(metrics.f1_score(y, y_dummyPred,average=None)).mean()
        # print("Dummy, most frequent acc:",dummy_freq_acc)

        # dummy_stratifiedRandom = DummyClassifier(strategy='stratified',random_state=0)
        # dummy_strat2= '{:.3%}'.format(metrics.accuracy_score(y, Get_yPred(X,y,clf_class=dummy_frequent))) #,sample_weight=balance_weights(y)))
        # 'print("Dummy, Stratified Random:",dummy_strat2)'
        print()

        resDict['dummy_freq:Accuracy'][fileName]=dummy_freq_acc
##        resDict['dummy_freq:f1'][fileName]=dummy_freq_f1 dummy_freq_f1_mean
        resDict['dummy_freq:f1'][fileName]=dummy_freq_f1_mean

        resDict['dummy_freq_f1_weighted'][fileName]=dummy_freq_f1_weighted
        # resDict.dummy_Stratfreq[fileName]=dummy_strat2

        "We can get seperately the best model for Acc, and the best for f1!"
        "WARNING!? In binary case - default F1 works for the 1 class, in sklearn 15. and lower"
        # bestEst_f1,bestScore_f1 = ModelParam_GridSearch(X,y,cv=3,scoreParam = 'f1')
        "Temporary workaround until next SKlearn update of F1 metric:"
        # bestEst_f1,bestScore_f1 = ModelParam_GridSearch(X,y,cv=3,scoreParam = 'f1')f1_scorer
        bestEst_f1,bestScore_f1 = ModelParam_GridSearch(X,y,cv=3,scoreParam = f1_scorer)

        bestEst_acc,bestScore_acc = ModelParam_GridSearch(X,y,cv=2,scoreParam = 'accuracy')
        print("bestEst (f1):",bestEst_f1)#,"best f1",bestScore_f1)
        print("bestEst (f1):",bestEst_acc)#,"best acc",bestScore_acc)

        #Temp
        # bestEst_f1=bestEst_acc=bestEst = RandomForestClassifier(n_jobs=-1)

        if GetRFEPerf==True:
            bestEst_RFE,bestScore_RFE = ModelParam_GridSearch(X_RFE,y,cv=3,scoreParam = 'f1')

        "Modified to get 2 estimators"
        scores_acc = cross_val_score(estimator=bestEst_acc, X=X, y=y, cv=StratifiedShuffleSplit(y, n_iter=13, test_size=0.18), n_jobs=-1) #Accuracy
        print("Accuracy: %0.3f (+- %0.2f)" % (scores_acc.mean(), scores_acc.std() * 2))
        scores_f1 = cross_val_score(estimator=bestEst_f1, X=X, y=y, cv=StratifiedShuffleSplit(y, n_iter=13, test_size=0.18), n_jobs=-1, scoring='f1')
        print("f1: %0.3f (+- %0.2f)" % (scores_f1.mean(), scores_f1.std() * 2))

        resDict['Accuracy'][fileName]=round(scores_acc.mean(),4)
        resDict['Accuracy_SD'][fileName]=round(scores_acc.std(),4)
        resDict['f1'][fileName]=round(scores_f1.mean(),4)
        resDict['f1_SD'][fileName]=round(scores_f1.std(),4)
        resDict['Array-f1-Scores'][fileName]=(scores_f1)
        resDict['Array-Acc-Scores'][fileName]=(scores_acc)
        resDict['bestML-f1'][fileName]=(str(bestEst_f1))
        resDict['bestML-Acc'][fileName]=(str(bestEst_acc))

        #ORIG
        # Acc,Acc_SD,f1,f1_SD = CV_multi_stats(X, y, bestEst,n=15)

        # resDict['Accuracy'][fileName]=round(Acc,4)
        # resDict['Accuracy_SD'][fileName]=round(Acc_SD,4)
        # resDict['f1 score'][fileName]=round(f1,4)
        # resDict['f1_SD'][fileName]=round(f1_SD,4)
        # resDict['Best (f1) Model parameters'][fileName]= bestEst

        print()
        # print(fileName," Done")

    print("Saving results to file")
    resDict.to_csv("OutputData.tsv", sep=',')


# #alt
# def plot_confusion_matrix(cm,target_names, title='Confusion matrix', cmap=plt.cm.Blues):
#     '''
#     http://scikit-learn.org/dev/auto_examples/model_selection/plot_confusion_matrix.html
#     NOT As Good!
#     '''

#     plt.imshow(cm, interpolation='nearest', cmap=cmap)
#     plt.title(title)
#     plt.colorbar()
#     tick_marks = np.arange(len(target_names))
#     plt.xticks(tick_marks, target_names, rotation=45)
#     plt.yticks(tick_marks, target_names)
#     plt.tight_layout()
#     plt.ylabel('True label')
#     plt.xlabel('Predicted label')

##############################################################################

if __name__ == '__main__' :

    GetAllPerf ()


# ################################################################################

# ##    os.chdir(r'E:\Dropbox\Dropbox\bioInf_lab\AA_info\fastas\Benchmarks/Thermophiles/')
#     # os.chdir(r'E:\Dropbox\Dropbox\bioInf_lab\AA_info\CODE\feat_extract\test_seq\NP\SP_CleavedNP+Neg')

#     # os.chdir(r'/a/fr-05/vol/protein/danofer/ProtFeat/feat_extract/test_seq/Thermophile')
#     os.chdir(r'/cs/prt3/danofer/ProtFeat/feat_extract/test_seq/NP/SP_Cleaved+NP+Neg_Big')

#     df = pd.read_csv('trainingSetFeatures.csv')
# ##    df.drop('proteinname',axis=1, inplace=True)
#     feature_cols = [col for col in df.columns if col not in ['classname','Id','proteinname']]
#     X=df[feature_cols].values
#     y=df.classname.values

#     le = LabelEncoder()
#     y = le.fit_transform(y)

#     Fwe = SelectFwe(alpha=0.01).fit(X,y)
#     X=Fwe.transform(X)
#     # feature_cols=feature_cols[Fwe.get_support()]
#     print("F filter: ",X.shape)


#     ModelParam_GridSearch(X, y, cv=3,scoreParam = 'roc_auc')

# #     from functools import partial
# #     from sklearn.metrics import precision_score, make_scorer
# #     import numpy as np
# # ##    custom_scorer = make_scorer(partial(precision_score, pos_label=1))
# #     custom_scorer = make_scorer(partial(f1_score, pos_label=1))

# #     clf,bestScore = ModelParam_GridSearch(X,y,cv=3,scoreParam = custom_scorer)

# #     cm = confusion_matrix(y, Get_yPred (X,y,clf))
# #     np.set_printoptions(precision=2)
# #     target_names=list(le.classes_)
# #     PlotConfusionMatrix (cm,target_names)

# #     # http://scikit-learn.org/dev/auto_examples/model_selection/plot_confusion_matrix.html
# #     cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
# #     print('Normalized confusion matrix')
# #     print(cm_normalized)
# #     plt.figure()
# #     plot_confusion_matrix(cm_normalized,target_names=list(le.classes_), title='Normalized confusion matrix')

# ################################################################################
##############################################################################

   # y_pred = Get_yPred(X,y,clf_class=bestEst,n_folds=4)
        # print("matthews_corrcoef",matthews_corrcoef(y,y_pred))
        # print("classification_report \n",classification_report(y,y_pred,target_names=lb_encoder.classes_,sample_weight=balance_weights(y)))
        # cm = metrics.confusion_matrix(y, y_pred)
        # print("confusion matrix: \n",cm)



    "Print out nice confusion matrices for classifiers:"
    "http://nbviewer.ipython.org/github/rasbt/musicmood/blob/master/code/classify_lyrics/nb_init_model.ipynb#Confusion-matrix"
    "  http://bugra.github.io/work/notes/2014-11-22/an-introduction-to-supervised-learning-scikit-learn/"

    'Decision tree visualization: http://scikit-learn.sourceforge.net/dev/modules/generated/sklearn.tree.export_graphviz.html#sklearn.tree.export_graphviz'
    'http://nbviewer.ipython.org/gist/greglandrum/4316460'

    'http://nbviewer.ipython.org/github/gmonce/scikit-learn-book/blob/master/Chapter%204%20-%20Advanced%20Features%20-%20Feature%20Engineering%20and%20Selection.ipynb'
