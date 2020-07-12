from feature_selection_utils import *
from joblib import Parallel, delayed

def find_best_params(X, y, params_grid, record_fi):

    if os.path.exists(record_fi):
        os.remove(record_fi)

    acc_record = pd.DataFrame(columns=["tree_num","max_depth","min_samples_split","acc",  "sensitivity", "specificity" ] )
    with open(record_fi, "a+") as fout:
        all_params = list(ParameterGrid(params_grid) )
        #all_params_df = pd.DataFrame(ParameterGrid(params_grid) )
        res = Parallel(n_jobs=24)(delayed(loo_validation)(params_['tree_num'],
                                                   params_["max_depth"],
                                                   params_['min_samples_split'], X, y) for params_ in all_params)
        acc_all, train_acc_all, sensitivity_all, specificity_all = zip(*res)
        for acc, train_acc, sensitivity, specificity,params in zip(acc_all, train_acc_all, sensitivity_all, specificity_all,all_params):
            #acc, train_acc, sensitivity, specificity = loo_validation(params["tree_num"], params["max_depth"],
            #                                                          params["min_samples_split"], X, y)
            line = "In tree_num={:d}, max_depth={:d}, min_samples_split={:d}, train_acc={:.4f}, loo_acc={:.4f}, sensitivity={:.4f}, specificity={:.4f}".format(
                params["tree_num"], params["max_depth"], params["min_samples_split"], train_acc, acc, sensitivity,
                specificity)
            #print(line)
            fout.write(line + "\n")

            acc_record_ = pd.DataFrame({"tree_num": [params["tree_num"]],
                                        "max_depth": [params["max_depth"]],
                                        "min_samples_split":[params["min_samples_split"]],
                                        "acc": [acc],
                                        "sensitivity": [sensitivity],
                                        "specificity": [specificity] } )
            acc_record= acc_record.append(acc_record_, ignore_index=True)

        '''
        acc_record_2 = pd.DataFrame({"tree_num": all_params_df["tree_num"],
                                        "max_depth": all_params_df["max_depth"],
                                        "min_samples_split": all_params_df["min_samples_split"],
                                        "acc": acc_all,
                                        "sensitivity": sensitivity_all,
                                        "specificity": specificity_all } )
        '''

    acc_record = acc_record.sort_values(by=[ 'acc'], ascending=False)
    acc_record = acc_record.reset_index(drop=True)
    #acc_record_2 = acc_record_2.sort_values(by=["sensitivity", "specificity", 'acc'], ascending=False, ignore_index=True)

    best_tree_num = acc_record.loc[0,"tree_num"]
    best_max_depth = acc_record.loc[0,"max_depth"]
    best_min_samples_split = acc_record.loc[0,"min_samples_split"]
    best_acc = acc_record.loc[0,"acc"]

    acc_record.to_csv("rf_my_record.csv", index=False)
    #acc_record_2.to_csv("rf_my_record_2.csv", index=False)
    print(best_tree_num, best_max_depth, best_min_samples_split, best_acc)
    return best_tree_num, best_max_depth, best_min_samples_split, best_acc


def save_feature_importance(X, y, feature_name, tree_num, max_depth, min_samples_split, fout):
    feature_importance = get_feature_importance(tree_num, max_depth, min_samples_split, X, y)
    feature_importance_dict = {"feature_name": feature_name, "feature_importance": feature_importance}
    feature_importance_df = pd.DataFrame.from_dict(feature_importance_dict)
    feature_importance_df.to_csv(fout, index=None)
    return feature_importance


tdir = "./FeatureSelection/"
chck_dir(tdir)
data_assemble = ["top50", "top100", "top200", "proc"]
#data_assemble  = ["proc"]
#data_assemble = ["top50", "top100", "top200"]

data_dir = ["./data/%s/" % s for s in data_assemble]
record_fi = [tdir + "%s_RF_params_setting_record.txt" % [s, "all"][s == "proc"] for s in data_assemble]
feature_imp_fi = [tdir + "feature(%s)_importance.csv" % [s, "all"][s == "proc"] for s in data_assemble]

params_grid = {"tree_num": [3, 5, 7, 10, 20, 30, 50, 100],
               "max_depth": [2, 3, 4, 7, 10],
               "min_samples_split": [2, 3, 4, 5, 7]
               }


best_params_record = {"data_dir": [], "tree_num": [],
                      "max_depth": [], "min_samples_split": [],
                      "best_loo_acc": []}

topk_for_eval = [10, 20, 50, 100, 200]

topk_for_eval_record = []

for i in range(len(data_assemble)):

    X, y, feature_name = load_data(data_dir[i])
    res = find_best_params(X, y, params_grid, record_fi[i] )
    # res = [10,3,4,1]
    best_tree_num = res[0]
    best_max_depth = res[1]
    best_min_samples_split = res[2]
    best_acc = res[3]

    best_params_record["data_dir"].append([data_assemble[i], "all"][data_assemble[i] == "proc"])
    best_params_record["tree_num"].append(best_tree_num)
    best_params_record["max_depth"].append(best_max_depth)
    best_params_record["min_samples_split"].append(best_min_samples_split)
    best_params_record["best_loo_acc"].append(best_acc)

    feature_imp = save_feature_importance(X, y, feature_name,
                            best_tree_num, best_max_depth,
                            best_min_samples_split, feature_imp_fi[i])

    # eval classification ability of top k features computed by RF
    # feature_imp = load_feature_imp(feature_imp_fi[i])
    for topk in topk_for_eval:
        X_imp, new_feature_name = select_top_imp_feature(X, feature_imp, feature_name, topk)
        fout = tdir + "feature(%s)_importance_top%d_for_retraining.csv" \
               % ([data_assemble[i], "all"][data_assemble[i] == "proc"], topk)
        acc, train_acc,sensitivity, specificity = eval_topk_RF_feature(X_imp, y, new_feature_name, best_tree_num,
                                              best_max_depth, best_min_samples_split, fout)

        line = "In {:s}, set tree_num={:d}, max_depth={:d}, min_samples_split={:d}, " \
               "select top{:d}, train_acc={:.4f}, loo_acc={:.4f}, sensitivity={:.4f}, specificity={:.4f}".format(
                feature_imp_fi[i], best_tree_num, best_max_depth, best_min_samples_split,
                topk, train_acc, acc,sensitivity, specificity)

        topk_for_eval_record.append(line)

best_params_record_df = pd.DataFrame.from_dict(best_params_record)
best_params_record_df.to_csv(tdir+"RF_best_params_settings_for_feature_selection.csv", index=False)

topk_for_eval_fout = tdir + "eval_RF_topk_features.txt"

if os.path.exists(topk_for_eval_fout):
    os.remove(topk_for_eval_fout)

topk_for_eval_record = [line + "\n" for line in topk_for_eval_record]
with open(topk_for_eval_fout, "a+") as f:
    f.writelines(topk_for_eval_record)


