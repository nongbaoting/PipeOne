import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, ParameterGrid
from sklearn.metrics import accuracy_score
import os
from sklearn.metrics import confusion_matrix

def chck_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir,exist_ok=True)

def sensitivity_specificity(y_true, y_pred):
    """
    计算敏感度和特异度
    :param y_true: 样本真实标签，(sample_num, ) array
    :param y_pred: 预测的样本标签, (sample_num, ) array
    :return: sensitivity, specificity
    """
    if sum(y_true) == len(y_true) and sum(y_pred) == len(y_pred):
        sensitivity = 1
        specificity = 0
    elif sum(y_true) == 0 and sum(y_pred) == 0:
        sensitivity = 0
        specificity = 1
    else:
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        sensitivity = 0
        specificity = 0
        if (tp+fn) > 0:
            sensitivity = tp / (tp + fn)
        if (tn+fp) > 0:
            specificity = tn / (tn + fp)

    return sensitivity, specificity


def loo_validation(tree_num, max_depth, min_samples_split, X, y):
    loo = LeaveOneOut()
    count = 0
    y_predict_lst = []
    y_test_lst = []
    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        clf = RandomForestClassifier(n_estimators=tree_num,
                                     max_depth=max_depth,
                                     min_samples_split=min_samples_split, random_state=0, n_jobs=24)
        clf.fit(X_train, y_train)
        y_predict = clf.predict(X_test)
        if 1.0 * y_predict[0] == 1.0 * y_test[0]:
            count += 1
        y_predict_lst.append(y_predict[0])
        y_test_lst.append(y_test[0])

    y_test_arr = np.array(y_test_lst)
    y_predict_arr = np.array(y_predict_lst)

    clf_train = RandomForestClassifier(n_estimators=tree_num,
                                       max_depth=max_depth,
                                       min_samples_split=min_samples_split, random_state=0, n_jobs= 24)
    clf_train.fit(X, y)
    y_train_pred = clf_train.predict(X)
    train_acc = accuracy_score(y, y_train_pred)

    #test_sensitivity, test_specificity = sensitivity_specificity(y_test_arr, y_predict_arr)
    test_sensitivity, test_specificity =0, 0
    acc = count / (1.0 * X.shape[0])
    return acc, train_acc, test_sensitivity, test_specificity


def get_feature_importance(tree_num, max_depth, min_samples_split, X, y):
    clf = RandomForestClassifier(n_estimators=tree_num,
                                 max_depth=max_depth,
                                 min_samples_split=min_samples_split, random_state=0)
    clf.fit(X, y)
    return clf.feature_importances_


def load_data(data_dir):
    data_fi = os.listdir(data_dir)
    X_lst = []
    feature_name = []
    y = None
    for fi in data_fi:
        if "clustering" in fi:
            y_info = pd.read_csv(data_dir + fi, header=0, index_col=None)
            y = np.array(y_info["Subtype"])
            #
            #y = np.array(y_info['subtype'].replace({'Primary Tumor':1, 'Solid Tissue Normal':0 }))

        else:
            data = pd.read_csv(data_dir+fi, header=0, index_col=0)
            datavalue = data.values
            feature_name = feature_name + list(data.columns)
            X_lst.append(datavalue)

    X_tuple = tuple(X_lst)
    X = np.column_stack(X_tuple)

    return X, y, feature_name


def load_feature_imp(fi):
    feature_imp_df = pd.read_csv(fi, index_col=None, header=0)
    feature_imp = np.array(feature_imp_df["feature_importance"].values)

    return feature_imp


def select_top_imp_feature(X, feature_imp, feature_name, topk):
    sort_idx = np.argsort(-feature_imp)
    nonzero_num = np.sum(feature_imp != 0)
    topk = np.minimum(topk, nonzero_num)
    X_imp = X[:, sort_idx[0:topk]]
    new_feature_name = [feature_name[i] for i in list(sort_idx[0:topk])]
    return X_imp, new_feature_name


def eval_topk_RF_feature(X, y, feature_name, tree_num, max_depth, min_samples_split, fout):

    # X_imp, new_feature_name = select_top_imp_feature(X, feature_imp, feature_name, topk)
    acc, train_acc,sensitivity, specificity = loo_validation(tree_num, max_depth, min_samples_split, X, y)
    feature_importance = get_feature_importance(tree_num, max_depth, min_samples_split, X, y)

    feature_importance_dict = {"feature_name": feature_name, "feature_importance": feature_importance}
    feature_importance_df = pd.DataFrame.from_dict(feature_importance_dict)
    # feature_importance_df.to_csv("feature(top100)_importance_top50_for_retraining.csv", index=None)
    feature_importance_df.to_csv(fout, index=None)

    return acc, train_acc,sensitivity, specificity
    # line = "In {:s}, set tree_num={:d}, max_depth={:d}, min_samples_split={:d}, select top{:d}, train_acc={:.4f}, " \
    #        "loo_acc={:.4f}".format(feature_imp_fi, tree_num, max_depth,
    #                                min_samples_split, topk, train_acc, acc)
    # print(line)
    # return line
