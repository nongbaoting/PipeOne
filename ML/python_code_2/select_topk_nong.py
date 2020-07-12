#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/6/26 11:11
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire,os,re
import numpy as np
import pandas as pd


def get_topk(X, weight, topk):
    idx = np.argsort(-weight)
    topk_idx = idx[:topk]
    return X.iloc[:, topk_idx]

def get_matrix_files(data_dir):
    view = []
    for entry in os.scandir(data_dir):
        if re.search('.csv$', entry.name):
            the_view = re.sub('^proc_', '', entry.name)
            the_view = re.sub('\.csv$', '', the_view)
            view.append(the_view)
    return view

class MYTOP:

    def select_topK(self, topK='50,100,200', consistence_file = "record_nmi.csv"):
        consist = pd.read_csv(consistence_file)
        consist = consist.sort_values(by='consistency', ascending=False )
        low_dim, alpha, gamma = consist.iloc[0,0:3]
        topk_ls = [int(i) for i in topK.split(',') ]
        proc_file_lst = get_matrix_files("./data/proc/")

        cluster =  "./results/lowDim=%d/" \
                    "lowDim=%d_alpha=%.2f_gamma=%.2f_clustering.csv" % (low_dim, low_dim, alpha, gamma)
        for topk in topk_ls:

            if not os.path.exists("./data/top%d/" % topk):
                os.makedirs("./data/top%d/" % topk)
            os.system( f'cp {cluster} ./data/top{topk}/')
            for file in proc_file_lst:
                fi = "./data/proc/proc_%s.csv" % file
                X = pd.read_csv(fi, header=0, index_col=0)

                weight_fi = "./results/lowDim=%d/weight_%s_lowDim=%d_alpha=%.2f_gamma=%.2f.csv" \
                            % (low_dim, file, low_dim, alpha, gamma)
                weight = pd.read_csv(weight_fi, header=None, index_col=None)
                weight_arr = np.array(weight[1])
                topk_X = get_topk(X, weight_arr, topk)

                fout = "./data/top%d/top%d_proc_%s.csv" % (topk, topk, file)
                topk_X.to_csv(fout, index=True)

if __name__ == '__main__':
    fire.Fire(MYTOP)