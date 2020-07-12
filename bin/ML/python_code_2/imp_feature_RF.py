from main import loo_validation, load_data, get_feature_importance
import os
import numpy as np
import pandas as pd


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


data_dir = "./data/top100/"
feature_imp_fi = "feature(top100)_importance.csv"
tree_num = 100
max_depth = 3
min_samples_split = 4
topk = 50
X, y, feature_name = load_data(data_dir)
feature_imp = load_feature_imp(feature_imp_fi)
X_imp, new_feature_name = select_top_imp_feature(X, feature_imp, feature_name, topk)


acc, train_acc = loo_validation(tree_num, max_depth, min_samples_split, X_imp, y)
new_feature_importance = get_feature_importance(tree_num, max_depth, min_samples_split, X_imp, y)

feature_importance_dict = {"feature_name": new_feature_name, "feature_importance": new_feature_importance}
feature_importance_df = pd.DataFrame.from_dict(feature_importance_dict)
feature_importance_df.to_csv("feature(top100)_importance_top50_for_retraining.csv", index=None)

line = "In {:s}, set tree_num={:d}, max_depth={:d}, min_samples_split={:d}, select top{:d}, train_acc={:.4f}, " \
       "loo_acc={:.4f}".format(feature_imp_fi, tree_num, max_depth,
                           min_samples_split, topk, train_acc, acc)
print(line)

# record_fi = "top_k_features_for_retraining.txt"
# with open(record_fi, "a+") as fout:
#     fout.write(line+"\n")

