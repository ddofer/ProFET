"""
=======================================
Receiver Operating Characteristic (ROC)
=======================================

Example of Receiver Operating Characteristic (ROC) metric to evaluate
classifier output quality.

ROC curves typically feature true positive rate on the Y axis, and false
positive rate on the X axis. This means that the top left corner of the plot is
the "ideal" point - a false positive rate of zero, and a true positive rate of
one. This is not very realistic, but it does mean that a larger area under the
curve (AUC) is usually better.

The "steepness" of ROC curves is also important, since it is ideal to maximize
the true positive rate while minimizing the false positive rate.

ROC curves are typically used in binary classification to study the output of
a classifier. In order to extend ROC curve and ROC area to multi-class
or multi-label classification, it is necessary to binarize the output. One ROC
curve can be drawn per label, but one can also draw a ROC curve by considering
each element of the label indicator matrix as a binary prediction
(micro-averaging).

.. note::

    See also :func:`sklearn.metrics.roc_auc_score`,
             :ref:`example_plot_roc_crossval.py`.

"""
import sys
# sys.path += [r'C:\Anaconda\Lib', r'C:\Anaconda\Lib\site-packages']
sys.path += [r'/cs/prt3/danofer/anaconda3',r'/cs/prt3/danofer/anaconda3/Lib',r'/cs/prt3/danofer/anaconda3/Lib/site-packages']

print(__doc__)
from sklearn.metrics import roc_curve, auc
from sklearn import metrics
from sys import argv
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression,  RandomizedLogisticRegression
from sklearn.svm import SVC,LinearSVC
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import StratifiedKFold,cross_val_score,StratifiedShuffleSplit

from sklearn.feature_selection import RFE, RFECV, SelectFdr,SelectFwe,SelectKBest
from sklearn.linear_model import RandomizedLogisticRegression

from sklearn.metrics import precision_recall_fscore_support,matthews_corrcoef, classification_report,confusion_matrix

import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.metrics import f1_score
from PipeTasks import Get_yPred,balance_weights
from collections import Counter
from sklearn.dummy import DummyClassifier
import os, fnmatch
import seaborn as sns
from sklearn.preprocessing import LabelEncoder

import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
import os, fnmatch
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
from PipeTasks import Get_yPred,balance_weights

# Import some data to play with
#########################################
os.chdir(r'/a/fr-05/vol/protein/danofer/ProtFeat/feat_extract/test_seq/Thermophile')
##os.chdir(r'/cs/prt3/danofer/ProtFeat/feat_extract/test_seq/NP/SP_Cleaved+NP+Neg_Big')

df = pd.read_csv('trainingSetFeatures.csv')
##    df.drop('proteinname',axis=1, inplace=True)
feature_cols = [col for col in df.columns if col not in ['classname','Id','proteinname']]
X=df[feature_cols].values
y=df.classname.values

Fwe = SelectFwe(alpha=0.01).fit(X,y)
X=Fwe.transform(X)

le = LabelEncoder()
y = le.fit_transform(y)
# Binarize the output
# y = label_binarize(y, classes=[0, 1, 2])
# y = label_binarize(y)

##n_classes = y.shape[1]
n_classes=len(set(y))
target_names=list(le.classes_)
print ("n_classes",n_classes,"target_names",target_names)
# shuffle and split training and test sets
##X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,
##                                                    random_state=0)

'Best hyper param  classifiers For Thermophiles:'
# svm=SVC(C=50, cache_size=1900, class_weight='auto', gamma=0.0, kernel='rbf', probability=True)
rf = RandomForestClassifier(n_jobs=-1,min_samples_split= 2, max_features= 0.4, min_samples_leaf= 2, n_estimators= 20, max_depth= 8)

'Best hyper param  classifiers For NP:'
##svm =SVC(C=50, cache_size=1500, class_weight='auto',gamma=0.0, kernel='rbf', max_iter=-1, probability=True )
##rf = RandomForestClassifier(max_depth= None, min_samples_split= 3, min_samples_leaf= 1, n_estimators= 250,  n_jobs= -1, max_features= auto)


# Learn to predict each class against the other
##classifier = OneVsRestClassifier(svm.SVC(kernel='rbf', probability=True,
##                                  cache_size=1200,random_state=random_state))
##y_score = classifier.fit(X_train, y_train).decision_function(X_test) #ORIG
#Use Y_pred here, maybe with dec_func?

# y_score = Get_yPred (X,y,svm,n_folds=4, pred_proba=True)

y_score = Get_yPred (X,y,rf,n_folds=2, pred_proba=True)

# Compute ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    # fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i]) #ORIG
    fpr[i], tpr[i], _ = roc_curve(y_true =y, y_score=y_score)
    roc_auc[i] = auc(fpr[i], tpr[i])

# Compute micro-average ROC curve and ROC area
# fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
fpr["micro"], tpr["micro"], _ = roc_curve(y.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
'''
# Plot of a ROC curve for a specific class
plt.figure()
plt.plot(fpr[2], tpr[2], label='ROC curve (area = %0.2f)' % roc_auc[2])
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.show()
'''
# Plot ROC curve
plt.figure()
plt.plot(fpr["micro"], tpr["micro"],
         label='micro-average ROC curve (area = {0:0.2f})'
               ''.format(roc_auc["micro"]))
for i in range(n_classes):
    # plt.plot(fpr[i], tpr[i], label='ROC curve of class {0} (area = {1:0.2f})'
    #                                ''.format(i, roc_auc[i]))
    plt.plot(fpr[i], tpr[i], label='ROC curve of class {0} (area = {1:0.2f})'.format(i, target_names[i]))


plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('multi-class Receiver operating characteristic')
plt.legend(loc="lower right")
sns.despine()
plt.ion()

plt.show()
