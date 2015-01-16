from __future__ import division
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
from sklearn.linear_model import LogisticRegression, SGDClassifier,RandomizedLogisticRegression

from numba.decorators import jit
from time import time
from scipy.stats import randint as sp_randint
from sklearn.grid_search import RandomizedSearchCV,GridSearchCV
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import SelectPercentile, f_classif, chi2, \
    SelectFwe, SelectFdr

from operator import itemgetter
from sklearn.pipeline import Pipeline

import pandas as pd

SEED = 4  # always use a seed for randomized procedures

N_TRIALS = 5
'''
BUG: Overlapping classnames lacks handling.
multiindex file handling not wokring currently. (workaround: manually delete first column, and drop label 'label' name)
    UserWarning: The sum of true positives and false positives are equal to zero for some labels.
    Precision is ill defined for those labels [0 1 2 3].
    The precision and recall are equal to zero for some labels.
    fbeta_score is ill defined for those labels [0 1 2 3].
    average=average)
'''

#==============================================================================
' http://nbviewer.ipython.org/github/tdhopper/Research-Triangle-Analysts--Intro-to-scikit-learn/blob/master/Intro%20to%20Scikit-Learn.ipynb '

"Get Imp. Feature's NAMES:"
'http://stackoverflow.com/questions/22361781/how-does-sklearn-random-forest-index-feature-importances?rq=1'
#==============================================================================

mean_auc = 0.0
mean_precision = 0.0
mean_recall = 0.0
mean_accuracy = 0

def reset_means() :
    global mean_auc, mean_precision, mean_recall, mean_accuracy
    mean_auc = 0.0
    mean_precision = 0.0
    mean_recall = 0.0
    mean_accuracy = 0


# 'Old. loads to numpy arrays'
@jit
def import_npArray_TrainingData() :
    fname_NEG = r"E:\Dropbox\Dropbox\Books & Resources\NPID Lab BioInformatics\NPID Lab Stuff\NPID Additional Datasets\Training\Standard_13\NEG\Features_NEG-.txt"
    fileObject_POS = open(fname_POS, 'r')
    test_POS = np.loadtxt(fileObject_POS, delimiter="\t")
    fileObject_NEG = open(fname_NEG, 'r')
    test_NEG = np.loadtxt(fileObject_NEG, delimiter="\t")
    #For label /classification/ generation: (Classify NP+ as 1, NP- as 0..):
    samples_POS_training = len(test_POS)
    samples_NEG_training = len(test_NEG)
    y_labels = ([1] * samples_POS_training + [0] * samples_NEG_training)
    y_labels = np.asarray(y_labels)  #make it a numpy array (needed for scikit )
    trainingSets = np.vstack((test_POS, test_NEG))
    return trainingSets, y_labels

# @jit
def load_data(filename="Feat_normalized.csv") :
    '''
    Load training data from csv file.  Load labels from it.
    Return matrix, training labels, encoder for labels.
    http://blog.yhathq.com/posts/predicting-customer-churn-with-sklearn.html
    http://stackoverflow.com/questions/21589177/using-multiple-features-with-scikit-learn?rq=1

    Labels could just be the names? : http://stackoverflow.com/questions/13300160/non-integer-class-labels-scikit-learn?rq=1
    '''
    df = pd.read_csv(filename, index_col=0)
    lb = LabelEncoder()
    labels = lb.fit_transform((df.index.values))

    print ("labels: %s %s" %(type(labels),labels))
    features = df.values
    # labels = LabelEncoder.transform(np.asarray(df['labels'].values))
    'This could be done more elegantly. Check index num for later filtering!!'
    'TODO: Is pop needed? List of col.values??  '

    feature_names=df.columns.values  #No pop. (nd array, no labels index here)
    print("%s features: " % (len(feature_names)))

    # classes = label_encoder.transform(np.asarray(df['labels']))
    print('encoded labels: %s' % (set(labels)))
    # print("feature_names: %s" %(feature_names))
    return (features, labels, lb,feature_names)

"For case of multiindex .. (file/class name (in 'labels' column) + sample name (unamrked col."
def multilabels_load_data(filename="Feat_normalized.csv") :
    '''
    Load training data from csv file.  Load labels from it.
    Return matrix, training labels, encoder for labels.
    http://blog.yhathq.com/posts/predicting-customer-churn-with-sklearn.html
    http://stackoverflow.com/questions/21589177/using-multiple-features-with-scikit-learn?rq=1

    Labels could just be the names? : http://stackoverflow.com/questions/13300160/non-integer-class-labels-scikit-learn?rq=1
    '''
    df = pd.read_csv(filename, index_col=1) #Note Col!
    lb = LabelEncoder()

    # print('classes: %s' %(df['labels'].values))
    labels = lb.fit_transform((df.index.values))
    # labels = lb.fit_transform((df['labels'].values))
    df.drop('labels', inplace = True)
    # print(df.iloc[0])

    print ("labels: %s %s" %(type(labels),labels))
    features = df.values
    features=np.delete(features,0,1) #Remove index.
    # print(features)

    # labels = LabelEncoder.transform(np.asarray(df['labels'].values))
    'Check these next bits - getting feature names.. '
    'This could be done more elegantly. Check index num for later filtering!!'
    "feat_names = df.columns.drop('Unnamed: 0') " # MAy be better - check ordering!
    feature_names = list(df.columns) # list may not be needed

    # feature_names=df.columns.values  #No pop. (nd array)
    # feature_names.pop(0)  #Remove index. #ORIG
    print("%s features: " % (len(feature_names)))
    # data = np.array(df['feature1'])
    # classes = label_encoder.transform(np.asarray(df['labels']))
    print('encoded labels: %s' % (set(labels)))
    # print("feature_names: %s" %(feature_names))
    return (features, labels, lb,feature_names)

#     ###################################################
#     'TEST labels!'
#     df = pd.read_csv(filename, index_col=0)
#     labels = df.index.values

#     print("labels: %s , type %s" %(labels, type(labels)))
#     training_data = df.values
#     training_data=scale_data(training_data)
#     # feature_names = list(df.columns)
#     # s = SVC(class_weight="auto", cache_size=2200, shrinking=True)
#     s= LinearSVC(penalty='l2', loss='l2', dual=False, class_weight='auto')
#     # s.fit(training_data, labels)
#     print ("string labels fitted")
#     print("CV:")
#     scores_f1 = cross_validation.cross_val_score(estimator=s, X=training_data, y=labels, cv=2,
#                                                  n_jobs=1, scoring='f1')
#     print("Model f1: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
#     # CV_multi_stats(training_data,labels,s)
#     'End TEST'
# ###################################################

#    @jit
def CV_Binary_stats(X, y, model) :
    '''
    http://scikit-learn.org/stable/modules/model_evaluation.html#classification-metrics
    Note that some of the metrics here ONLY work for BINARY tasks.

    http://scikit-learn.org/stable/modules/cross_validation.html#computing-cross-validated-metrics
    '''
    global mean_auc, mean_precision, mean_recall, mean_accuracy
    n = N_TRIALS  # repeat the CV procedure 10 times to get more precise results
    for i in range(n) :
        # for each iteration, randomly hold out 30% of the data as CV set
        X_train, X_cv, y_train, y_cv = cross_validation.train_test_split(X, y,
                                                                         test_size=.3,
                                                                         random_state=i * SEED)
        # train model and make predictions
        model.fit(X_train, y_train)
        preds = model.predict(X_cv)
        fpr, tpr, thresholds = metrics.roc_curve(y_cv, preds)
        roc_auc = metrics.auc(fpr, tpr)
        print("( %d/%d)" % (i + 1, n))
        mean_auc += roc_auc
        accuracy = accuracy_score(y_cv, preds)
        precision = precision_score(y_cv, preds)
        recall = recall_score(y_cv, preds)
        mean_accuracy += accuracy
        mean_precision += precision
        mean_recall += recall
    mean_accuracy = (mean_accuracy / n)
    mean_precision = mean_precision / n
    mean_recall = mean_recall / n
    mean_auc = mean_auc / n
    print('mean_accuracy:  %s ' %(round(mean_accuracy, 2)))
    print('mean_precision:  %s ' %(round(mean_precision, 2)))
    print('mean_recall:  %s ' %(round(mean_recall, 2)))
    print('mean_auc:  %s ' %(round(mean_auc, 2)))
    reset_means()


def CV_multi_stats(X, y, model) :
    '''
    http://scikit-learn.org/stable/modules/model_evaluation.html#classification-metrics
    This version uses multiclass (or multilabel) compatible metrics.

    May be expanded to use the cross_val_score helper function:
    http://scikit-learn.org/stable/modules/generated/sklearn.cross_validation.cross_val_score.html
    http://scikit-learn.org/stable/modules/cross_validation.html#computing-cross-validated-metrics
    '''
    n = N_TRIALS  # repeat the CV procedure 10 times to get more precise results
    scores = cross_validation.cross_val_score(estimator=model, X=X, y=y) #, cv=n) #Accuracy
    scores_f1 = cross_validation.cross_val_score(estimator=model, X=X, y=y, scoring='f1')
    print("Model Accuracy: %0.3f (+- %0.2f)" % (scores.mean(), scores.std() * 2))
    print("Model f1: %0.3f (+- %0.2f)" % (scores_f1.mean(), scores_f1.std() * 2))


'TODO: Implement.  (More plots and test figures in source)'
def Feature_Importance_plot (est,names):
    'http://nbviewer.ipython.org/github/pprett/pydata-gbrt-tutorial/blob/master/gbrt-tutorial.ipynb'
    fx_imp = pd.Series(est.feature_importances_, index=names)
    fx_imp /= fx_imp.max()  # normalize
    fx_imp.sort()
    fx_imp.plot(kind='barh', figsize=FIGSIZE)

# @jit
def report(grid_scores, n_top=3) :
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores) :
        print("Model with rank: {0}".format(i + 1))
        print("Mean validation score: {0:.3f} (std: {1:.3f})".format(score.mean_validation_score,
                                                                     np.std(score.cv_validation_scores)))
        print("Parameters: {0}".format(score.parameters))
        print("")


from scipy import stats
from scipy.stats import randint as sp_randint
Tree_param_dist = {"n_estimators" : sp_randint(50, 500),
                   "max_depth" : sp_randint(2, 13), "max_features" : sp_randint(5, 70),
                   "min_samples_split" : sp_randint(1, 8), "min_samples_leaf" : sp_randint(1, 8),
                   "criterion" : ["gini"], "bootstrap" : [True]}
# "criterion": ["gini", "entropy"]}

SVM_param_grid = [{'C' : [0.1,0.001, 1, 10, 100, 1000], 'kernel' : ['linear', 'rbf'], 'gamma' : [0.001, 0.0001]}]
SVM_param_dist = [{'C' : stats.expon(scale=100), 'gamma' : stats.expon(scale=.1),
                   'kernel' : ['rbf', 'linear'], 'class_weight' : ['auto']}]

SGD_param_dist = [{'n_iter':stats.expon(scale=10) ,'learning_rate': ['optimal', 'invscaling'],
                   'l1_ratio' : stats.expon(scale=.1), 'penalty':['l2', 'l1' , 'elasticnet'],
                   'class_weight' : ['auto'] }]
SGD_param_grid = [{'n_iter':[10, 70,550, 2000,180000] ,'learning_rate': ['optimal', 'invscaling'],
   'l1_ratio' : [0.15,0.3,0.6], 'penalty':['l1' ,'elasticnet'],
   'loss': ["hinge","modified_huber","log"],'class_weight':['auto'] }]

def GridParamSearch(param_dist, clf, X, y, n_iter_search=15) :
    '''
    Searches using rand.search for best model paramters
    diff paramters searched by model type..
    http://nbviewer.ipython.org/github/treycausey/thespread/blob/master/notebooks/basic_random_forest_wp_model.ipynb?create=1
    @param clf: estimator/predictor used.
    @param param_dist: Grid of Parameter ranges to tune for the predictor,
    using randomized CV search.
    '''
    print("Starting grid parameter search")
    random_search = RandomizedSearchCV(clf, param_distributions=param_dist,
                                       n_iter=n_iter_search,n_jobs=-1)
    start = time()
    # random_search.fit(features, target)
    random_search.fit(X, y)
    print("RandomizedSearchCV took %.2f seconds for %d candidates"
          " parameter settings." % ((time() - start), n_iter_search))
    report(random_search.grid_scores_)


def scale_data(X) :
    scale_scaler = MinMaxScaler()
    normal_scaler = StandardScaler()
    # X=scaler.fit(X).transform(X)
    X = normal_scaler.fit_transform(X)
    X = scale_scaler.fit_transform(X)

    return X

'Ignore this bit - just a mess of calls to CV methods. '
def call_GridParamSearch_featfilt(X, y) :
    '''
        (def is Currently just a cut & paste from "main".)
        Calles def GridParamSearch , (which uses randomized CV to find odel param)
    Used to try different ml models, then get their optimal paramters
    '''
    print("SPARSE (L1) EXT gridparam scores:")
    #   clf = Pipeline([
    #       ('feature_selection', LinearSVC(penalty="l1", loss='l1',dual=False, class_weight='auto')),
    # ('classification', ExtraTreesClassifier(n_jobs=3)
    #   )])
    'Sparse; L1 penalized features selection prior to RF fitting/prediction'
    clf_svm = LinearSVC(penalty="l1", loss='l2', dual=False, class_weight='auto')
    clf_logit = LogisticRegression(penalty="l1", dual=False, class_weight='auto')

    'http://scikit-learn.org/0.13/auto_examples/plot_feature_selection.html'
    print('Original features matrix:')
    print(X.shape)
    # Univariate feature selection with F-test for feature scoring
    # We use the default selection function: the 20% most significant features
    # selector = SelectPercentile(f_classif, percentile=20)
    selector = SelectPercentile(chi2, percentile=20)
    X_anova = selector.fit_transform(X, y)
    print(
        'New (2 f_classif) Using statistical feature selection: features matrix is:')
    print(X_anova.shape)

    # lda = LDA(n_components=10)
    # X_lda = lda.fit_transform(X, y)
    # print('New LDA filtered features matrix:')
    # print(X_lda.shape)

    X_svm = clf_svm.fit_transform(X, y)  #Get Sparse feature selections..
    # print(clf.feature_importances_ )
    print('New sparse (SVM filtered) features matrix:')
    print(X_svm.shape)

    print("Res of SVM fitting of (F scores filtered =2) for more feature selection:")
    X_doubleFilt_svm_f = clf_svm.fit_transform(X_anova, y)
    print(X_doubleFilt_svm_f.shape)
    print("param search on sparse features matrix")
    GridParamSearch(param_dist=Tree_param_dist, clf=clf_EXT, X=X_svm, y=y)

def plot_BestKFeatures (X_train, y_train):
    '''
    http://nbviewer.ipython.org/github/gmonce/scikit-learn-book/blob/master/Chapter%204%20-%20Advanced%20Features%20-%20Feature%20Engineering%20and%20Selection.ipynb
    Find the best percentile of features to use,
    using cross-validation on the training set and get K best feats
    '''
    from sklearn import cross_validation
    from sklearn import feature_selection
    from sklearn import tree
    dt = tree.DecisionTreeClassifier(criterion='entropy')
    dt = RandomForestClassifier(n_jobs=2, bootstrap=True, n_estimators=250, criterion='gini')
    dt = dt.fit(X_train, y_train)

    percentiles = range(1, 95, 5)
    results = []
    for i in range(1, 95, 5):
        fs = feature_selection.SelectPercentile(feature_selection.chi2, percentile=i) #Original
        fs = feature_selection.SelectPercentile(feature_selection.f_classif, percentile=i) # alt
        X_train_fs = fs.fit_transform(X_train, y_train)
        scores = cross_validation.cross_val_score(dt, X_train_fs, y_train, cv=4)
        #print i,scores.mean()
        results = np.append(results, scores.mean())

    optimal_percentil = np.where(results == results.max())[0]
    print (("Optimal number of features:{0}".format(percentiles[optimal_percentil])), "\n")

    # Plot number of features VS. cross-validation scores
    import pylab as pl
    import matplotlib.pylab as pl
    pl.figure()
    pl.xlabel("Number of features selected")
    pl.ylabel("Cross validation accuracy)")
    pl.plot(percentiles,results)
    print ("Mean scores:",results)
    return

def plot_RFE(X,y):
    from sklearn.svm import SVC
    from sklearn.cross_validation import StratifiedKFold
    from sklearn.feature_selection import RFECV
    from sklearn.datasets import make_classification
    from sklearn.metrics import zero_one_loss
    import pylab as pl
    import matplotlib.pylab as pl

    # Create the RFE object and compute a cross-validated score.
    # svc= SVC(kernel="linear", class_weight="auto", cache_size=1200, shrinking=True)
    svc=LinearSVC(penalty='l1', loss='l2', dual=False, class_weight='auto',multi_class='ovr')
#    SGD = SGDClassifier(penalty='elasticnet',class_weight='auto',n_jobs=-1,n_iter=10,l1_ratio =0.15)
##    rfecv = RFECV(estimator=svc, step=0.1, cv=StratifiedKFold(y, 5), scoring='roc_auc')
    rfecv = RFECV(estimator=svc, step=0.2,cv=StratifiedKFold(y, 2), scoring='f1')
    X_RFE = rfecv.fit_transform(X, y)

    print("Optimal number of features in X_RFE : %d" % rfecv.n_features_)
    # Plot number of features VS. cross-validation scores
    pl.figure()
    pl.xlabel("Number of features selected")
    pl.ylabel("Cross validation score (nb of misclassifications)")
    pl.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
    pl.show()
    print ('RFE Opt.shapes features CV score:')
    CV_multi_stats(X_RFE,y,svc)
    return (X_RFE,rfecv)


if __name__ == '__main__' :
    # filename = r".\training_data\mammal_locations\mod_Organelle_Mammal_Loc_normalized.csv"
    # filename = r".\training_data\Extracellular\ExtraCell_Train_Feat.csv"
    # filename = r".\training_data\Viri-Bact-Vert\Vir-0.9\Feat_norm.csv"
    # filename = r".\training_data\Viri-Capsids\Feat.csv"
#    filename = r".\training_data\NP\Feat.csv"
#    filename = r".\training_data\NP\NP+SP_Small\Norm_Feat.csv"
    filename = r".\training_data\hsp33\Norm_Feat.csv"

    features, labels, lb_encoder,feature_names = load_data(filename)
    # features, labels, lb_encoder = multilabels_load_data(filename)

    X, y = features, labels
    class_names=lb_encoder.inverse_transform(y)
    print("Data and labels imported")
    print(X.shape)
    # print('X: %s' %(X))
    # print('y: %s' %(y))

#    'Normalizing Unneeded if features already normalized. (scaling may still be needed)'
    X = scale_data(X)
    print("Features Data scaled")

#    SGD = SGDClassifier(penalty='elasticnet',class_weight='auto',n_jobs=-1,n_iter=35,l1_ratio =0.2)
    svc = LinearSVC(class_weight='auto')
    model_rf = RandomForestClassifier(n_jobs=-1, bootstrap=True, n_estimators=180,
                                        min_samples_leaf=3, min_samples_split =3,
                                        criterion='gini',compute_importances=True, max_depth=6)

    SVC_RBF= SVC(kernel="rbf", class_weight="auto", cache_size=2600, shrinking=True)
    SVC_linear= SVC(kernel="poly", cache_size=2700, shrinking=True)


    # model_rf.fit(X,y)
    # X_SGD = model_rf.transform(X, threshold='1.5*mean') # forests!
    X_SGD = model_rf.fit_transform(X,y)
    print('X Reduced (by RF) features amount:')
    print(X_SGD.shape)

    def ReducedFeaturesDF(X,y):
        '''
        Returns a dataframe with only a subset of features/columns retained
        '''
        from sklearn.feature_selection import RFE
        est = LinearSVC( penalty='l1', loss='l2', dual=False, class_weight='auto')
#        selectK = SelectKBest(score_func = f_classif, k=45)
        selectRFE = RFE(estimator=est, n_features_to_select=22, step=0.15)
        selectK=selectRFE

        selectK.fit(X,y)
        selectK_mask=selectK.get_support()
        K_featnames = feature_names[selectK_mask]
        print ("reduced RFE features:")
        print(K_featnames)
        Reduced_df = pd.read_csv(filename, index_col=0)
        Reduced_df = Reduced_df[Reduced_df.columns[selectK_mask]]
#        Reduced_df.to_csv('REDUCED_Feat.csv')
        return Reduced_df

#    ReducedFeaturesDF(X,y)
    # z=pd.DataFrame(data=X_SGD,index=y)
    # z.to_csv('REDUCED_Feat.csv')

    '! TODO: Are feature names being saved consistently??? Check with get K best'
    print('Plotting best features percent. 1)X_SGD ')
##    plot_BestKFeatures (X_SGD, y)
##    print('Plotting best features percent. 2) Full X')
##    plot_BestKFeatures (X, y)

#    est = LinearSVC(penalty='l1', loss='l2', dual=False, class_weight='auto')
#    selectRFE = RFE(estimator=est, n_features_to_select=30, step=0.1)
#    selectRFE.fit(X,y)
    print("Reduce feat DF")

    selectK = SelectKBest(k=25)
    selectK.fit(X,y)
    selectK_mask=selectK.get_support()
    K_featnames = feature_names[selectK_mask]
    print("K_featnames: %s" %(K_featnames))
    k2=SelectKBest(k=60)
    X_bestK=k2.fit_transform(X_SGD,y)

    # print("runA SGD:")
    # selectK = SelectKBest()
    # selectK.fit(X_SGD,y)
    # selectK_mask=selectK.get_support()
    # K_featnames = feature_names[selectK_mask]
    # print("K_featnames SGD: %s" %(K_featnames))

    '''
    X_RFE,rfecv = plot_RFE(X,y)
    rfe_mask = rfecv.get_support() # Boolean mask.
    rfe_featnames = feature_names[rfe_mask]
    print("RFE-cv selected features from X (orig.) are: %s" %(rfe_featnames))
    '''

    "following Probably doesn't work!! Due to num columns changing.."
#    X_SGD_RFE,SGD_rfecv = plot_RFE(X_SGD,y)
#    rfe_SGD_featnames = feature_names[SGD_rfecv.get_support()]
#    print("RFE selected features from X-SGD (reduced) are: %s" %(rfe_SGD_featnames))

##    plot_RFE(X,y)
    # print('RFE Plotted')
    'At this stage we can also, easily, get a new DF from the data. (fun for plotting, saving..).'
    'df_RFE = pd.DataFrame(X_RFE, index=classes, columns=rfe_featnames)'

    '''Another option for getting feature names + importances: '
    'From: http://www.datarobot.com/blog/python-getting-started-with-data-science/   '
    xtrain, xtest, ytrain, ytest = cross_validation.train_test_split( X,y, train_size=0.4)
    importances = pd.Series(gbm.feature_importances_, index=data.columns)
    print importances.order(ascending=False)[:10]

    'We can also plot with radviz on a new df, (once we reset the index):'
    df_RFE = pd.DataFrame(X_RFE, index=classes, columns=rfe_featnames)
    df_RFE.index.rename('class',inplace=True)
    df_RFE.reset_index(inplace=True)
    import pylab
    radviz(df_RFE,'class')
    pylab.show()
    '''

    print('Getting CV stats (For initially reduced matrix):')

    print('X CV stats: - predicting with RF Classifier')
    CV_multi_stats(X,y,model_rf)
    print('X_SGD CV stats: - predicting with RF Classifier')
    CV_multi_stats(X_SGD,y,model_rf)
    # print('X_bestK CV stats: - predicting with RF Classifier')
    # CV_multi_stats(X_bestK,y,model_rf)


    # rbf_model = SVC(C=10, gamma=0.5, kernel="rbf", class_weight="auto", cache_size=3200, shrinking=True)
    # linear_model= SVC(C=10, gamma=0.5, kernel="rbf", class_weight="auto", cache_size=3200, shrinking=True)
    linear_model = LinearSVC( penalty='l2', loss='l2', dual=False, class_weight='auto')

    'http://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_regularization.html '
    clf_GBT = GradientBoostingClassifier(n_estimators=1500, learning_rate=0.05,
                                         max_depth=4,subsample=0.45,min_samples_leaf=5)

    # print("ExtraTreesClassifier gridparam scores:")
    # GridParamSearch(param_dist=Tree_param_dist,clf=clf_EXT,X=X,y=y)


    # print("classes: %s" %(list(lb_encoder.classes_)))
    # print ('Random Forest stats:')
    # CV_multi_stats(X,y,model_rf)
    # # # CV_Binary_stats (X,y,model_rf)

    clf_SVC = SVC(class_weight="auto", cache_size=3200, shrinking=True)

    # selector_fdr = SelectFdr(alpha=0.01)
    # X_fdr = selector_fdr.fit_transform(X, y)
    # clf_logit = LogisticRegression(penalty="l1", dual=False, class_weight='auto')
    # clf_svm_L2 = LinearSVC(penalty="l2", loss='l1', dual=True, class_weight='auto')


    # print('Logit reduced features matrix:')
    # print(X_new_logit.shape)


