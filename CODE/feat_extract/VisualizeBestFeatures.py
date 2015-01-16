
# coding: utf-8

# In[1]:

from sklearn.preprocessing import LabelEncoder
import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import StratifiedKFold,cross_val_score,StratifiedShuffleSplit
from sklearn.feature_selection import RFE,SelectKBest,SelectFwe
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.preprocessing import balance_weights
import numpy as np
from sklearn import feature_selection, tree  # "http://nbviewer.ipython.org/github/gmonce/scikit-learn-book/blob/master/Chapter%204%20-%20Advanced%20Features%20-%20Feature%20Engineering%20and%20Selection.ipynb"
pd.set_option('display.max_columns', 15)
pd.set_option('display.max_rows', 5)
# mpl.rc('title', labelsize=6)
mpl.rc('ytick', labelsize=6)
mpl.rc('xtick', labelsize=3)

class RandomForestClassifierWithCoef(RandomForestClassifier):
    '''
    Let's us use RF with RFE feature selection.
    http://stackoverflow.com/questions/24123498/recursive-feature-elimination-on-random-forest-using-scikit-learn?rq=1
    '''

    def fit(self, *args, **kwargs):
        super(RandomForestClassifierWithCoef, self).fit(*args, **kwargs)
        self.coef_ = self.feature_importances_


if __name__ == '__main__' :

    # In[2]:

    # os.chdir(r'E:\Dropbox\Dropbox\bioInf_lab\AA_info\fastas\Benchmarks/Thermophiles/')
    os.chdir(r'/a/fr-05/vol/protein/danofer/ProtFeat/feat_extract/test_seq/Thermophile')
    dataName = 'Thermophiles'

    # os.chdir(r'E:\Dropbox\Dropbox\bioInf_lab\AA_info\CODE\feat_extract\test_seq\NP\SP_CleavedNP+Neg')
    # os.chdir(r'/cs/prt3/danofer/ProtFeat/feat_extract/test_seq/NP/SP_Cleaved+NP+Neg_Big')
    # dataName = 'Neuropeptide Precursors'

    # In[3]:

    df = pd.read_csv('trainingSetFeatures.csv')
    # df.drop('proteinname',axis=1, inplace=True)

    # In[ ]:

    feature_cols = [col for col in df.columns if col not in ['classname','Id','proteinname']]
    feature_cols=np.array(feature_cols)


    # In[ ]:

    X=df[feature_cols].values
    y=df.classname.values

    # In[ ]:
    le = LabelEncoder()
    y = le.fit_transform(y)


    # In[ ]:
    Fwe = SelectFwe(alpha=0.01).fit(X,y)
    X=Fwe.transform(X)
    print("F-test -> ",X.shape)
    feature_cols=feature_cols[Fwe.get_support()]
    # feature_cols
    # In[]:


    # In[ ]:

    from operator import itemgetter
    def report(grid_scores, n_top=1) :
        '''
        Print out top models/parameters after a grid search for model params.
        '''
        top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
        for i, score in enumerate(top_scores) :
            if n_top>1:
                print("Model with rank: {0}".format(i + 1))
            print("Average Cross-Validation score (while tuning): {0:.3f} (std: {1:.2f})".format(
                score.mean_validation_score, np.std(score.cv_validation_scores)))
            print("Model Parameters: {0}".format(score.parameters))
            print("")


    # # In[ ]:

    param_dist = {"max_depth": [6,8, None],
                  "max_features": ['auto',0.4],
                  "min_samples_leaf": [1,2],
                  "bootstrap": [True, False],
                  'min_samples_split':[2,3],
                  "criterion": [ "gini"],
                  "n_estimators":[200],
                  "n_jobs":[-1]}
    # rf = RandomForestClassifierWithCoef(n_estimators=150,n_jobs=-1)
    rf = RandomForestClassifierWithCoef(max_depth= 7, min_samples_split= 3, min_samples_leaf= 2, n_estimators= 160,  n_jobs= -1, max_features= "auto")
    # gs = GridSearchCV(rf, param_grid=param_dist,
    #     n_jobs=-2,cv=3, iid=False,scoring='precision',
    #     fit_params={'sample_weight': balance_weights(y)})
    # gs.fit(X, y)

    # # # In[ ]:

    # report(gs.grid_scores_, n_top=2)
    # # #Use best param forest but with more estimators - to do feature pruning.
    # rf = gs.best_estimator_

    # In[8]:
    # rf = RandomForestClassifierWithCoef(max_depth= None, min_samples_split= 3, min_samples_leaf= 1, n_estimators= 150,  n_jobs= -1, max_features= "auto")


    # # In[ ]:


    scores = cross_val_score(rf,X,y,n_jobs=-2,cv=StratifiedShuffleSplit(y,n_iter=11,test_size=0.2))
    print("X RF Accuracy: %0.3f (+- %0.2f)" % (scores.mean(), scores.std() * 2))
    scores_f1 = cross_val_score(rf,X,y,n_jobs=-2,cv=StratifiedShuffleSplit(y,n_iter=11,test_size=0.2),scoring='f1')
    print("X RF f1: %0.3f (+- %0.2f)" % (scores_f1.mean(), scores_f1.std() * 2))


    # # In[11]:

    # rf.fit(X,y,sample_weight= balance_weights(y))


    # # # In[12]:

    # x2 = rf.transform(X,'7*mean')
    # print(x2.shape)


    # # In[ ]:

    # x2_scores=(cross_val_score(rf,x2,y,n_jobs=-1,cv=StratifiedShuffleSplit(y,n_iter=8,test_size=0.25),scoring='f1')).mean()
    # print("x2_scores:",x2_scores)
    # print("With just",x2.shape[1],"features, we have %f performance!" %(100*x2_scores/scores_f1.mean()))

    # # In[ ]:

    # from sklearn.svm import LinearSVC
    # svc = LinearSVC(C=10, penalty='l1', dual=False)
    # svc.fit(X, y)
    # selected_feature_names = feature_cols[[list(set(np.where(svc.coef_ != 0)[-1]))]]

    # In[ ]:

    rfeSelect = RFE(estimator=rf,n_features_to_select=15, step=0.02)
    X_RFE = rfeSelect.fit_transform(X,y)
    print(X_RFE.shape)

    # In[ ]:

    RFE_FeatureNames = feature_cols[rfeSelect.get_support()]
    print(RFE_FeatureNames)


    # In[ ]:

    # gbc = GradientBoostingClassifier(n_estimators=100)
    # print("X_RFE - GB score:",cross_val_score(gbc,X_RFE,y,n_jobs=-1,cv=7).mean())
    # print("X - GB score:",cross_val_score(gbc,X,y,n_jobs=-1,cv=7).mean())

    # In[ ]:

    RFE_ScoreRatio = 100*(cross_val_score(rf,X_RFE,y,n_jobs=-1,cv=StratifiedShuffleSplit(y,n_iter=11,test_size=0.2),scoring='f1').mean())/scores_f1.mean()
    print("Even with just",X_RFE.shape[1]," features, we have %f performance!" %(RFE_ScoreRatio))


    # In[ ]:

    def PlotFeaturesImportance(X,y,featureNames):
        '''
        Plot the relative contribution/importance of the features.
        Best to reduce to top X features first - for interpretability
        Code example from:
        http://bugra.github.io/work/notes/2014-11-22/an-introduction-to-supervised-learning-scikit-learn/
        '''

        gbc = GradientBoostingClassifier(n_estimators=60)
        gbc.fit(X, y)
        # Get Feature Importance from the classifier
        feature_importance = gbc.feature_importances_
        # Normalize The Features
        feature_importance = 100 * (feature_importance / feature_importance.max())
        sorted_idx = np.argsort(feature_importance)
        pos = np.arange(sorted_idx.shape[0]) + 3.5
        # pos = np.arange(sorted_idx.shape[0])
        # plt.figure(figsize=(16, 12))
        plt.figure(figsize=(16, 12), dpi=250)
        plt.barh(pos, feature_importance[sorted_idx], align='center', color='#7A68A6')
        #plt.yticks(pos, np.asanyarray(df.columns.tolist())[sorted_idx]) #ORIG
        plt.yticks(pos, np.asanyarray(featureNames)[sorted_idx])

        plt.xlabel('Relative Importance')
        plt.title('%s: Top Features' %(dataName))
        plt.grid('off')
        plt.ion()
        plt.show()
        plt.savefig(str(dataName)+'TopFeatures.png',dpi=220)


    def altPlotFeaturesImportance(X,y,featureNames):
    "http://nbviewer.ipython.org/github/cs109/2014/blob/master/homework-solutions/HW5-solutions.ipynb"
        clf = RandomForestClassifier(n_estimators=60)

        clf.fit(X,Y)
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
    # In[ ]:



    PlotFeaturesImportance(X_RFE,y,RFE_FeatureNames)


    # In[ ]:




    # from sklearn.externals.six import StringIO
    # from sklearn import tree
    # clf = tree.DecisionTreeClassifier()



