
import numpy as np
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.metrics import mutual_info_score
import os


def subtype_to_digit(subtype):
    if subtype == "I":
        subtype_digit = 0
    elif subtype == "II":
        subtype_digit = 1
    else:
        subtype_digit = 2

    return subtype_digit


def eval(latentX, cluster_num, defined_subtype, sample_id, event, time_to_event, save_path):
    predict_subtype = KMeans(n_clusters=cluster_num,random_state=0).fit_predict(latentX)
    #nmi = mutual_info_score(defined_subtype, predict_subtype)
    nmi = 0
    cluster_dict = {"Run": sample_id, "event": event, "time_to_event": time_to_event,
                    "subtype": predict_subtype}
    cluster_df = pd.DataFrame.from_dict(cluster_dict)
    cluster_df.to_csv(save_path, header=True, index=False)
    return nmi


clinical_fi = "./data/sample.cli.csv"

clinical = pd.read_csv(clinical_fi, header=0, index_col=None)

'''
expr_subtype = list(clinical["expression_based_subtype"])
expr_subtype = [subtype_to_digit(subtype) for subtype in expr_subtype]
'''

sample_id = list(clinical["Run"])

event = clinical["event"].copy()
time_to_event = clinical["time_to_event"]
#event[event == "Last follow-up"] = 0
#event[event == "Disease progression"] = 1


event = list(event)
time_to_event = list(time_to_event)

alpha = [0.01, 0.02, 0.1, 1, 5]
gamma = [0, 0.1, 1, 10, 100]
low_dim = [2, 3, 4, 5, 6, 7]
# low_dim = [7]
cluster_num = 3

record_fout = "record_nmi.txt"

if os.path.exists(record_fout):
    os.remove(record_fout)

with open(record_fout, "a+") as fwrite:
    for low_dim_ in low_dim:
        for alpha_ in alpha:
            for gamma_ in gamma:
                latent_fi = "./results/lowDim=%d/" \
                            "X_lowDim=%d_alpha=%.2f_gamma=%.2f.csv" \
                            % (low_dim_, low_dim_, alpha_, gamma_)
                latentX_df = pd.read_csv(latent_fi, header=None, index_col=0)


                latentX = latentX_df.values

                save_path = "./results/lowDim=%d/" \
                            "lowDim=%d_alpha=%.2f_gamma=%.2f_clustering.csv" \
                            % (low_dim_, low_dim_, alpha_, gamma_)
                expr_subtype =1 ## 临时
                nmi = eval(latentX, cluster_num, expr_subtype, sample_id,
                           event, time_to_event, save_path)

                line = "In low_dim=%d, alpha=%.2f, gamma=%.2f, consistency between " \
                       "cluster subtype and expression based subtype: %.3f" \
                       % (low_dim_, alpha_, gamma_, nmi)
                fwrite.write(line+"\n")
                print(line)







