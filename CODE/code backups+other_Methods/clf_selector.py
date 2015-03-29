from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier
from sklearn.cross_validation import ShuffleSplit
import numpy as np
from sklearn.metrics import f1_score, roc_auc_score

def main():
    xtrain=np.load('data/x_train.npy')
    ytrain=np.load('data/y_train.npy')
    
    #train-test split
    ss=ShuffleSplit(np.shape(ytrain)[0],n_iter=1,test_size=0.3, random_state=42)
    for train_idx, test_idx in ss:
        xtest=xtrain[test_idx,:]
        ytest=ytrain[test_idx]
        xtrain=xtrain[train_idx,:]
        ytrain=ytrain[train_idx]
    
    clf_et=ExtraTreesClassifier(n_estimators=300,random_state=42,verbose=2)
    clf_et.fit(xtrain,ytrain)
    et_preds=clf_et.predict(xtest)
    print "initial f1-score, et-classifier:", f1_score(ytest,et_preds)
    feat_imp=clf_et.feature_importances_
    sorted_fi=feat_imp[np.argsort(feat_imp)[::-1]] #descending sort
    clf_gb=GradientBoostingClassifier(n_estimators=120, random_state=42)
    feats_tot=np.shape(xtrain)[1]
    f1_best=0
    print "output format:"
    print "no of features, f1-score, roc-score of class-predictions, roc-score of probabilities"
    for feats in range(1,feats_tot+1):
        threshold_idx=min(len(sorted_fi),feats)
        threshold=sorted_fi[threshold_idx]
        select=(feat_imp>threshold)
        clf_gb.fit(xtrain[:,select],ytrain)
        tmp_preds=clf_gb.predict(xtest[:,select])
        tmp_probs=clf_gb.predict_proba(xtest[:,select])[:,1]
        f1=f1_score(ytest,tmp_preds)
        roc_pred=roc_auc_score(ytest,tmp_preds)
        roc_prob=roc_auc_score(ytest,tmp_probs)
        if f1>f1_best:
            f1_best=f1
            np.save('features/clf_sel.npy',select)
        print feats,f1,roc_pred,roc_prob
        if feats>25:
            break
    print "f1_best:", f1_best
    
if __name__=="__main__":
    main()