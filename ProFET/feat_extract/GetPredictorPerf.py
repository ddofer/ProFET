
#cross_val_predict - taken from release 0.16.0
from sklearn.cross_validation import StratifiedKFold,cross_val_score,cross_val_predict
from sklearn import metrics
import pandas as pd
import numpy as np

# Use a constant seed
np.random.seed(1274)
SEED = 14  # always use a seed for randomized procedures



def get_scores(scores,y,label=None, verbose=True):
    '''
    Returns a dictionary of metrics for a given classification of the data (given by Cross_val_predict).
    scores: list
        Classifier predictions on data
    y: list
        True Class labels
    label: string
        Name of the classifier used
    '''
    results_dict = {}
    try:
        roc_auc_no_avg = metrics.roc_auc_score(y, scores,average=None)
        if verbose: print("roc_auc for each class: %0.4f " % (roc_auc_no_avg))
        results_dict['ROC_AUC not averaged'] = roc_auc_no_avg
    
        roc_auc_weighted = metrics.roc_auc_score(y, scores,average='weighted')
        if verbose: print("roc_auc (weighted-Av): %0.4f " % (roc_auc_weighted))
        results_dict['ROC_AUC weighted'] = roc_auc_weighted

    except ValueError as e:
        print(e)

    f1_pos = metrics.f1_score(y, scores,average='binary')
    if verbose: print("POS f1: %0.4f  " % (f1_pos))
    results_dict['F1'] = f1_pos
    
    av_PR = metrics.average_precision_score(y, scores) # corresponds to the area under the precision-recall curve
    if verbose: print("Av_precision (Prec-Recall AUC): %0.3f " % (av_PR))
    results_dict['Averaged precision'] = av_PR
    
    accuracy = metrics.accuracy_score(y, scores)
    if verbose: print("Accuracy: %0.3f " % (accuracy))
    results_dict['Accuracy'] = accuracy
    
    precision,recall,fscore,support = metrics.precision_recall_fscore_support(y, scores,average='binary')
    if verbose: print("Precision: %0.3f " % (precision))
    results_dict['Precision'] = precision
    
    if verbose: print("Recall: %0.3f " % (recall))
    results_dict['Recall'] = recall
    
    # if verbose: print("fscore(fBeta): %0.4f  [%s]" % (fscore, label))
    mcc = metrics.matthews_corrcoef(y, scores)
    if verbose: print("MCC: %0.3f " % (mcc))
    results_dict['Matthews correlation coefficient'] = mcc

    results_dict = {k:round(float(v),4) for k, v in results_dict.items()}
    return results_dict

def demo_getPerf(X,y,Classifier,Classifier_label):
    """
    Classifier_type: Sklearn model
        Type of classifier to use, and it's parameters
    Classifier_label: string
        Descriptive Name of the classifier. e.g "Forest"
    """
    results = {}

    scores = cross_val_predict(Classifier, X, y, cv=10, n_jobs=-1)
    results[label] = get_scores(scores,y,Classifier_label)

    res_df = pd.DataFrame(results)
    res_df.to_csv(outputFileName+"tsv", sep='\t')

if __name__ == '__main__':
    from ProFET import *
    input_file = 'ProFET_features.csv'
    features, labels, label_encoder, feature_names = load_data(input_file)
    inputmodel_file = 'test_model'
    model, label_encoder, scaler, model_feature_names = load_model_from_file(inputmodel_file)
    Classifier = model
    X = features
    y = labels
    Classifier_label = 'forest'
    results = {}
    scores = cross_val_predict(Classifier, X, y, cv=10, n_jobs=-1)
