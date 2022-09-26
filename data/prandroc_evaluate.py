def auroc(test_list,test_true):
    import numpy as np
    import sklearn.metrics
    from sklearn.metrics import roc_curve, auc
    y_true = np.array(test_true)
    y_pred = np.array(test_list)
    fpr, tpr, thresholds = roc_curve(y_true, y_pred, pos_label=1)#pos_label设置的是什么标签是正的
    auroc = auc(fpr,tpr)
    return auroc,fpr,tpr

def aupr(test_list,test_true):
    from sklearn.metrics import precision_recall_curve,auc
    import numpy as np
    import sklearn.metrics
    y_true = np.array(test_true)
    y_list = np.array(test_list)
    precision ,recall,thresholds = precision_recall_curve(y_true,y_list)
    aupr = auc(recall,precision)
    return aupr,recall,precision

