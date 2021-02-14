import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
import os


def load_data(data_dir):
    data_fi = os.listdir(data_dir)
    X_lst = []
    feature_name = []
    y = None
    for fi in data_fi:
        data_full = os.path.join(data_dir, fi)
        if "clustering" in fi:
            y_info = pd.read_csv(data_full, header=0, index_col=None)
            y = np.array(y_info["Group"])
            #y = np.array(y_info["Group"].replace({'Primary Tumor': 1, 'Solid Tissue Normal':0 }))
        else:
            data = pd.read_csv(data_full, header=0, index_col=0)
            datavalue = data.values
            feature_name = feature_name + list(data.columns)
            X_lst.append(datavalue)

    X_tuple = tuple(X_lst)
    X = np.column_stack(X_tuple)

    return X, y, feature_name

class FeatureValidation:
    def __init__(self, X_train, y_train, X_test, y_test):

        self.X_train = X_train
        self.y_train = y_train
        self.X_test = X_test
        self.y_test = y_test

    def select_topk_feature(self, feature_importance, feature_name, topk):
        X_train = self.X_train
        X_test = self.X_test

        sort_idx = np.argsort(-feature_importance)

        ## use non-zero feature importance only
        nonzero_num = np.sum(feature_importance != 0)
        topk = np.minimum(topk, nonzero_num )
        
        topk_X_train = X_train[:, sort_idx[0:topk]]
        topk_X_test = X_test[:, sort_idx[0:topk]]
        topK_feature_name =  [feature_name[i] for i in  list(sort_idx[0:topk])  ]

        self.topk_X_train = topk_X_train
        self.top_X_test = topk_X_test
        self.topK_feature_name = topK_feature_name
        # return topk_X_train, topk_X_train

    def do_grid_search(self, params, cv, n_jobs, max_features = None, random_state = 42):
        X_train = self.topk_X_train
        y_train = self.y_train
        classifier = RandomForestClassifier(max_features= max_features, random_state= random_state)
        clf = GridSearchCV(classifier, params, n_jobs=n_jobs, cv=cv)

        clf.fit(X_train, y_train)

        best_params = clf.best_params_
        best_score = clf.best_score_
        best_estimator = clf.best_estimator_

        ## add training evaluation
        y_pred = best_estimator.predict(X_train )
        train_sen, train_spec = self.sensitivity_specificity(y_train, y_pred)

        return best_estimator, best_params, best_score, train_sen, train_spec

    def sensitivity_specificity(self, y_true, y_pred):
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        senitivity = 0
        specificity = 0
        if (tp+fn) > 0:
            senitivity = tp / (tp + fn)
        if (tn+fp) > 0:
            specificity = tn / (tn + fp)

        return senitivity, specificity

    def evaluation(self, estimator):
        X_test = self.top_X_test
        y_test = self.y_test
        y_pred = estimator.predict(X_test)
        test_sen, test_spec = self.sensitivity_specificity(y_test, y_pred)

        return test_sen, test_spec

    def get_feature_importance(self, fout, params, max_features = None, random_state = 42):
        clf = RandomForestClassifier(n_estimators=params["n_estimators"], max_depth=params["max_depth"],
                                 min_samples_split=params["min_samples_split"],
                                 class_weight=params["class_weight"], max_features = max_features, random_state=random_state) 

        clf.fit(self.topk_X_train, self.y_train)
        feature_importance = clf.feature_importances_
        feature_dict = {"feature": self.topK_feature_name, "weight": feature_importance}

        feature_df = pd.DataFrame.from_dict(feature_dict)
        feature_df.to_csv(fout, header=True, index=False)
        #return feature_importance


def grid_search(X, y, classifier, parameters, cv, n_jobs):
    clf = GridSearchCV(classifier, parameters, n_jobs=n_jobs, cv=cv)
    clf.fit(X, y)
    best_params = clf.best_params_
    best_score = clf.best_score_
    return best_params, best_score


def get_feature_importance(X, y, feature_name, fout, params, random_state = 42):
    clf = RandomForestClassifier(n_estimators=params["n_estimators"], max_depth=params["max_depth"],
                                 min_samples_split=params["min_samples_split"],
                                 class_weight=params["class_weight"], random_state=random_state)

    clf.fit(X, y)
    feature_importance = clf.feature_importances_
    feature_dict = {"feature": feature_name, "weight": feature_importance}

    feature_df = pd.DataFrame.from_dict(feature_dict)
    feature_df.to_csv(fout, header=True, index=False)
    return feature_importance














