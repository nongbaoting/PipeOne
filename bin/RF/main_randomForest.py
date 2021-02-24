from feature_validation import *
from sklearn.model_selection import LeaveOneOut
import os,fire
from collections import defaultdict

def main(threads = 12, train_dir = "./data/train_dir/", test_dir ="./data/test_dir/", 
    n_estimators=[3,5,7,10,20,30,50,100], max_depth=[2,3,4,7,10], min_samples_split = [2,3,4,5,7],
    max_features = 'auto', random_state = 42):

    logfi = "./log.txt"  # log file, record key points in the whole process
    if os.path.exists(logfi):
        os.remove(logfi)

    f = open(logfi, "a+")

    data_dir = "./data/"  # where data files located, set accordingly.
    n_jobs = threads

    '''
    # Please set accordingly before running the code
    X_train = None  # features in training set, (sample_num, feature_num) array
    y_train = None  # binary numerical label in training set, e.g. 0 stands for tumor and 1 stands for normal, (sample_num, ) array
    feature_name = None # feature name for all features in training set, list, string, e.g. [feature_name1, feature_name2, ...]
    X_test = None # features in test set, (sample_num, feature_num) array
    y_test = None # the same as y_train, except it is for test set
    '''

    X_train, y_train, feature_name = load_data(train_dir)
    X_test , y_test , feature_name_t = load_data(test_dir )
    X_train, X_test = np.nan_to_num(X_train), np.nan_to_num(X_test)
    #print(y_train)
    print(n_estimators)
    parameters = {"n_estimators": [ int(i) for i in n_estimators ],
                "max_depth": [ int(i) for i in max_depth ],
                "min_samples_split": [ int(i) for i in min_samples_split ],
                "class_weight": [None, "balanced", "balanced_subsample"]
                }

    # select the best parameter settings for training set
    cv = LeaveOneOut()
    classifier = RandomForestClassifier(random_state= random_state)
    best_params, best_score = grid_search(X_train, y_train, classifier, parameters, cv, n_jobs)

    chck1 = "The best parameter settings found in grid search"

    f.write(chck1+"\n")
    print(chck1)
    print("best_params: " + str(best_params) )
    print("best_score: " + str(best_score))
    
    for key in best_params.keys():
        chck2 = key + "=" + str(best_params[key])
        f.write(chck2 + "\n")
        print(chck2)

    chck3 = "The best parameter settings achieved accuracy: %.5f" % best_score
    f.write(chck3+"\n")
    print(chck3)
    #
    # # save feature importance under the best parameter settings
    feature_importance_out = data_dir + "feature_importance.csv"
    feature_importance = get_feature_importance(X_train, y_train, feature_name, feature_importance_out, best_params, random_state = random_state)

    chck4 = "feature importance is saved in " + feature_importance_out
    f.write(chck4+"\n")


    feature_validation = FeatureValidation(X_train, y_train, X_test, y_test)

    topk_ls = [10, 20, 50, 100, 200, 500, len(feature_name) ]
    # topk_ls = [len(feature_name)]
    eval_dict = {"topk": topk_ls, "test_sensitivity": [], "test_specificity": [] }
    eval_dt = defaultdict(list)
    for topk in topk_ls:
        eval_dt["topK"].append(topk)
        # update topk features
        feature_validation.select_topk_feature(feature_importance, feature_name, topk)

        # find the best parameter settings with top k features in training set.
        #max_features = 'auto'
        val_best_estimator, val_best_params, val_best_score, val_train_sen, val_train_spec \
             = feature_validation.do_grid_search(parameters, cv, n_jobs, max_features =  max_features, random_state = random_state)
        
        val_feature_importance_out = data_dir + f"/feature_importance_{topk}.csv"
        feature_validation.get_feature_importance(val_feature_importance_out, val_best_params, 
                            max_features =  max_features, random_state = random_state)
        
        chck5 = "The best parameter settings with top%d features" % topk
        f.write(chck5+"\n")
        print(chck5)
        for key in val_best_params.keys():
            chck6 = key + "=" + str(val_best_params[key])
            f.write(chck6+"\n")
            print(chck6)
            eval_dt[key].append(str(val_best_params[key]))

        chck7 = "The best parameter settings achieved accuracy: %.5f" % val_best_score
        f.write(chck7+"\n")
        print(chck7)

        # evaluation with top k features in test set
        test_sensitivity, test_specificity = feature_validation.evaluation(val_best_estimator)

        eval_dict["test_sensitivity" ].append(test_sensitivity)
        eval_dict["test_specificity" ].append(test_specificity)

        eval_dt["train_sensitivity"].append(val_train_sen)
        eval_dt["train_specificity"].append(val_train_spec)
        eval_dt["test_sensitivity" ].append(test_sensitivity)
        eval_dt["test_specificity" ].append(test_specificity)

        # acc_ls.append(val_acc)
        # chck8 = "With top%d features, the validation set achieved accuracy:%.5f" % (topk, val_acc)
        # f.write(chck8+"\n")
        # print(chck8)

    # eval_df = pd.DataFrame.from_dict(eval_dict)
    # eval_df_fout = data_dir + "discriminative_power_of_topk_feature.csv"
    # eval_df.to_csv(eval_df_fout, index=False, header=True)

    eval_df2 = pd.DataFrame.from_dict( eval_dt )
    eval_df2_fout = data_dir + "discriminative_power_of_topk_feature.csv"
    eval_df2.to_csv(eval_df2_fout, index=False, header=True )

    f.write("Evalution metric is saved in " + eval_df2_fout + "\n")
    f.close()
    #
    # f.close()
    #
    # fig, ax = plt.subplots()
    # ax.plot(topk_ls, acc_ls)
    # ax.set(xlabel="Top K", ylabel="Accuracy in validation set",
    #        title="Discriminative power of top K features")
    # figout = "discriminative_power_of_top_k_features.png"
    # fig.savefig(figout)

if __name__ == "__main__":
    fire.Fire(main)








