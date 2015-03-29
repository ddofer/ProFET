# -*- coding: utf-8 -*-
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
from sklearn import feature_selection
from sklearn.svm import LinearSVC

class RandomForestClassifierWithCoef(RandomForestClassifier):
    '''
    Let's us use RF with RFE feature selection.
    http://stackoverflow.com/questions/24123498/recursive-feature-elimination-on-random-forest-using-scikit-learn?rq=1
    '''

    def fit(self, *args, **kwargs):
        super(RandomForestClassifierWithCoef, self).fit(*args, **kwargs)
        self.coef_ = self.feature_importances_


if __name__ == '__main__' :


#    os.chdir('/a/fr-05/vol/protein/danofer/ProtFeat/feat_extract/chap/train/')
    os.chdir('/cs/prt3/danofer/ProtFeat/feat_extract/test_seq/BenchM/scop/Scop_PartialClass_10')
    dataName = 'Scop-Partial_10'

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
    print("Orig X -> ",X.shape)
    Fwe = SelectFwe(alpha=0.001).fit(X,y)
    X=Fwe.transform(X)
    print("F-test -> ",X.shape)
    feature_cols=feature_cols[Fwe.get_support()]

# In[ ]:

    rf = RandomForestClassifierWithCoef(max_depth= 9, min_samples_split= 3, min_samples_leaf= 3, n_estimators= 650,  n_jobs= -1, max_features= "auto")


    # In[ ]:

    scores = cross_val_score(rf,X,y,n_jobs=-1,cv=StratifiedShuffleSplit(y,n_iter=7,test_size=0.3))
    print("X RF Accuracy: %0.3f (+- %0.2f)" % (scores.mean(), scores.std() * 2))
#    scores_f1 = cross_val_score(rf,X,y,n_jobs=-1,cv=StratifiedShuffleSplit(y,n_iter=10,test_size=0.22),scoring='f1')
#    print("X RF f1: %0.3f (+- %0.2f)" % (scores_f1.mean(), scores_f1.std() * 2))

     # In[ ]:
    svc = LinearSVC(C=20, penalty='l1', dual=False)
    svc.fit(X, y)
    selected_feature_names = feature_cols[[list(set(np.where(svc.coef_ != 0)[-1]))]]
    X_svm = svc.transform(X)
    print("X_svm L1 transformed:", X_svm.shape)
    X=X_svm


     # In[ ]:

    rfeSelect = RFE(estimator=rf,n_features_to_select=10, step=0.15)
    X_RFE = rfeSelect.fit_transform(X,y)
    print(X_RFE.shape)

    # In[ ]:

    RFE_FeatureNames = feature_cols[rfeSelect.get_support()]
    print("RFE_FeatureNames: \n",RFE_FeatureNames)


    # In[ ]:

    "http://stackoverflow.com/questions/21548750/plotting-histograms-against-classes-in-pandas-matplotlib"
    for featName in RFE_FeatureNames:
        df.groupby("class").feature.hist(alpha=0.4)
        df.groupby("classname")[featName].plot(kind='kde')


