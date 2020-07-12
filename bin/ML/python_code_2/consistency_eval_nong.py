#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/6/24 15:16
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'

import fire, os, re
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import mutual_info_score

def subtype_to_digit(subtype):
    if subtype == "I":
        subtype_digit = 0
    elif subtype == "II":
        subtype_digit = 1
    else:
        subtype_digit = 2

    return subtype_digit


def eval(latentX, cluster_num, defined_subtype, sample_id, save_path):
    predict_subtype = KMeans(n_clusters=cluster_num, random_state=0).fit_predict(latentX)
    nmi = mutual_info_score(defined_subtype, predict_subtype)
    cluster_dict = {"Run": sample_id, "defined_subtype": defined_subtype,
                    "subtype": predict_subtype}
    cluster_df = pd.DataFrame.from_dict(cluster_dict)
    cluster_df.to_csv(save_path, header=True, index=False)
    return nmi


def get_head_idx(df, pattern = re.compile("Run|Sample") ):
    for i, c in enumerate(df.columns):
        if pattern.match(c):
            return i

class MYCONSIST:

    def eval_2_groups(self, sample_group_file, cluster_num = 2):
        clinical = pd.read_csv(sample_group_file)
        sample_idx = get_head_idx(clinical, re.compile("Run|Sample") )
        group_idx  = get_head_idx(clinical, re.compile("Group") )

        sample_id = list( clinical.iloc[:, sample_idx] )
        subtype   = list( clinical.iloc[:, group_idx] )

        record_fout = "record_nmi.csv"

        if os.path.exists(record_fout):
            os.remove(record_fout)

        alpha = [0.01, 0.02, 0.1, 1, 5]
        gamma = [0, 0.1, 1, 10, 100 ]
        low_dim = [2, 3, 4, 5, 6, 7 ]

        params_res = pd.DataFrame(columns=['low_dim', 'alpha', 'gamma', 'consistency'])
        for low_dim_ in low_dim:
            for alpha_ in alpha:
                for gamma_ in gamma:
                    latent_fi = "./results/lowDim=%d/" \
                                "X_lowDim=%d_alpha=%.2f_gamma=%.2f.csv" \
                                % (low_dim_, low_dim_, alpha_, gamma_)
                    latentX_df = pd.read_csv(latent_fi, header=None, index_col=0)
                    latentX_df = latentX_df.reindex(sample_id)
                    latentX = latentX_df.values

                    save_path = "./results/lowDim=%d/" \
                                "lowDim=%d_alpha=%.2f_gamma=%.2f_clustering.csv" \
                                % (low_dim_, low_dim_, alpha_, gamma_)
                    nmi = eval(latentX, cluster_num, subtype, sample_id, save_path)
                    line = "In low_dim=%d, alpha=%.2f, gamma=%.2f, consistency between " \
                           "cluster subtype and expression based subtype: %.3f" \
                           % (low_dim_, alpha_, gamma_, nmi)
                    # fwrite.write(line + "\n")
                    print(line)

                    p_df = pd.DataFrame.from_dict({'low_dim': [low_dim_], 'alpha': [alpha_],
                                                   'gamma': [gamma_], 'consistency': [nmi] })
                    params_res = params_res.append(p_df, ignore_index=True )

        params_res.sort_values(by='consistency', ascending=False
                               ).to_csv("record_nmi.csv", header=True, index=False)


if __name__ == '__main__':
    fire.Fire(MYCONSIST)