#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/6/26 11:11
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire,os,re,sys 
import numpy as np
import pandas as pd

def chck_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def get_topk(X, weight, topk):
    idx = np.argsort(-weight)
    topk_idx = idx[:topk]
    return X.iloc[:, topk_idx]

def get_matrix_files(data_dir):
    view = []
    for entry in os.scandir(data_dir):
        if re.search('.csv$', entry.name) and not re.search('clustering.csv$', entry.name):
            the_view = re.sub('^proc_', '', entry.name)
            the_view = re.sub('\.csv$', '', the_view)
            view.append(the_view)
    return view

def split_params(NMF_param):
    low_dim_, alpha_, gamma_ = [i.split('=')[1] for i in NMF_param.split('_')]
    return (int(low_dim_), float(alpha_), float(gamma_) )

def split_fileName_params(fileName):
    low_dim, alpha, gamma, cluster_n = [i.split('=')[1] for i in fileName.split('_')[0:4]]
    return (int(low_dim), float(alpha),float(gamma),int(cluster_n) )

def select_topK( topK_importance='50,100,200', outdir = "./data_randomForest", cluster_survival_file = "record_log_rank_test_pvalue.csv", 
 cluster_file = ""):
    sv_fi = ''
    if os.path.exists(cluster_file):
        cluster_fi = cluster_file
        low_dim, alpha, gamma, cluster_n = split_fileName_params( os.path.basename(cluster_file) )
    elif os.path.exists(cluster_survival_file):
        consist = pd.read_csv(cluster_survival_file)
        consist = consist[consist.significant == True ].sort_values(by='mean_siloutte_width' , ascending=False)
        NMF_param, cluster = consist.iloc[0,0:2]
        low_dim, alpha, gamma = split_params(NMF_param)
        cluster_fi =  "./clusters/eval_cluster_num/%s_clusters=%d_clustering.csv" % (NMF_param, cluster)
        sv_fi      = "./clusters/surv_curve/%s_clusters=%d__clustering.pdf" % (NMF_param, cluster)
    else:
        sys.exit("please provide --cluster_survival_file or --cluster_file")

    chck_dir(outdir); 
    os.system(f'cp -r ./data/proc/  {outdir}/')
    os.system(f'cp {cluster_fi}  {outdir}/proc/')
    # copy selected cluster results
    os.system(f'cp {cluster_fi}  {outdir}/')
    if sv_fi:
        os.system(f'cp {sv_fi}  {outdir}/')
    proc_file_lst = get_matrix_files("./data/proc/")
    topk_ls = [int(i) for i in topK_importance.split(',') ]
    for topk in topk_ls:

        if not os.path.exists("%s/top%d/" % (outdir, topk) ) :
            os.makedirs("%s/top%d/" % (outdir, topk) )
        os.system( f'cp {cluster_fi} {outdir}/top{topk}/')
        for file_ in proc_file_lst:
            fi = "./data/proc/proc_%s.csv" % file_
            X = pd.read_csv(fi, header=0, index_col=0)

            weight_fi = "./NMF/lowDim=%d/weight_%s_lowDim=%d_alpha=%.2f_gamma=%.2f.csv" \
                        % (low_dim, file_, low_dim, alpha, gamma)
            weight = pd.read_csv(weight_fi, header=None, index_col=None)
            weight_arr = np.array(weight[1])
            topk_X = get_topk(X, weight_arr, topk)

            fout = "%s/top%d/top%d_proc_%s.csv" % (outdir, topk, topk, file_ )
            topk_X.to_csv(fout, index=True)

if __name__ == '__main__':
    fire.Fire(select_topK)