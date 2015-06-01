
from sklearn import metrics
from sys import argv
import numpy as np
import pandas as pd

from sklearn.cross_validation import StratifiedKFold,cross_val_score,StratifiedShuffleSplit
from operator import itemgetter
from sklearn.metrics import precision_recall_fscore_support,matthews_corrcoef, classification_report,confusion_matrix
from Model_trainer import load_data
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.metrics import f1_score
from PipeTasks import Get_yPred,balance_weights
from collections import Counter
from sklearn.dummy import DummyClassifier
import os, fnmatch
import seaborn as sns
from sklearn.preprocessing import LabelEncoder


def find_files(directory, pattern):
    'http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python'
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

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
    Dirr = r'K:\Lab-DataSets\test_seq'
    ##Dirr = r'./test_seq'
    if filePaths is None:
        filePaths = list(find_files(directory=Dirr, pattern='trainingSetFeatures.csv'))

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
         ,'bestML-Acc','bestML-f1'])


    #redDict holds results for each file/class, for saving to output-file

    i=-1
    for filePath in filePaths:
        i +=1

        'http://pythonconquerstheuniverse.wordpress.com/2008/06/04/gotcha-%E2%80%94-backslashes-in-windows-filenames/'
        filePath = os.path.normpath(filePath)
        print(filePath)
        fileName=str(fileNames[i])
        print()
        print("fileName: %s \n" %(fileName))
        "resDict['Name']= fileName"

        X, y, lb_encoder,featureNames = load_data(filePath, 'file') # X, y = features, labels

        print("Dummy classifiers output:")
        dummy_frequent = DummyClassifier(strategy='most_frequent',random_state=0)
        y_dummyPred = Get_yPred(X,y,clf_class=dummy_frequent)
        dummy_freq_acc = '{:.3}'.format(metrics.accuracy_score(y,y_dummyPred ))
        dummy_freq_f1_mean=(metrics.f1_score(y, y_dummyPred,average=None)).mean()
        dummy_freq_f1 = '{:.3}'.format(metrics.f1_score(y, y_dummyPred,average='weighted'))

        print ("f1 Av:",dummy_freq_f1_mean)
##        print ("f1 micro:",metrics.f1_score(y, y_dummyPred,average='micro'))
        print ("dummy_freq_f1",dummy_freq_f1)
        print("Dummy, most frequent acc:",dummy_freq_acc)
##        print("Diffence: dummy_freq_f1_def",dummy_freq_f1_def)
        print("\n classification_report: \n",metrics.classification_report(y, y_dummyPred))

        print()

        resDict['dummy_freq:Accuracy'][fileName]=dummy_freq_acc
        resDict['dummy_freq:f1'][fileName]=dummy_freq_f1_mean

    resDict.to_csv("dum.csv")

if __name__ == '__main__' :

    GetAllPerf ()

