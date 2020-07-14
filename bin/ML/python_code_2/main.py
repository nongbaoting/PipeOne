import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, ParameterGrid
from sklearn.metrics import accuracy_score
import os

def loo_validation(tree_num, max_depth, min_samples_split, X, y):
    loo = LeaveOneOut()
    count = 0
    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        clf = RandomForestClassifier(n_estimators=tree_num,
                               max_depth=max_depth,
                               min_samples_split=min_samples_split,random_state=0, n_jobs=4)
        clf.fit(X_train, y_train)
        y_predict = clf.predict(X_test)
        if 1.0*y_predict[0] == 1.0*y_test[0]:
            count += 1

    clf_train = RandomForestClassifier(n_estimators=tree_num,
                               max_depth=max_depth,
                               min_samples_split=min_samples_split,random_state=0, n_jobs=4)
    clf_train.fit(X, y)
    y_train_pred = clf_train.predict(X)
    train_acc = accuracy_score(y, y_train_pred)
    acc = count / (1.0*X.shape[0])
    return acc, train_acc

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
            y_info = pd.read_csv(data_dir+fi, header=0, index_col=None)
            y = np.array(y_info["subtype"])
        else:
            data = pd.read_csv(data_dir+fi, header=0, index_col=0)
            datavalue = data.values
            feature_name = feature_name + list(data.columns)
            X_lst.append(datavalue)

    X_tuple = tuple(X_lst)
    X = np.column_stack(X_tuple)

    return X, y, feature_name

data_dir = "./data/top200/"
tree_num = 30
max_depth = 4
min_samples_split = 7

X, y, feature_name = load_data(data_dir)

feature_importance = get_feature_importance(tree_num, max_depth, min_samples_split, X, y)

feature_importance_dict = {"feature_name": feature_name, "feature_importance": feature_importance}
feature_importance_df = pd.DataFrame.from_dict(feature_importance_dict)
feature_importance_df.to_csv("feature(top200)_importance.csv", index=None)
print()

# params_grid = {"tree_num": [3, 5, 7, 10, 20, 30, 50, 100],
#                "max_depth":[2, 3, 4, 7, 10],
#                "min_samples_split":[2, 3, 4, 5, 7]
#                }
#
# record_fi = "all_record.txt"
#
# with open(record_fi, "a+") as fout:
#     for params in ParameterGrid(params_grid):
#         # print(params)
#
#         acc, train_acc = loo_validation(params["tree_num"], params["max_depth"], params["min_samples_split"], X, y)
#         line = "In tree_num={:d}, max_depth={:d}, min_samples_split={:d}, train_acc={:.4f}, loo_acc={:.4f}".format(
#             params["tree_num"], params["max_depth"], params["min_samples_split"], train_acc, acc)
#         print(line)
#         fout.write(line+"\n")
