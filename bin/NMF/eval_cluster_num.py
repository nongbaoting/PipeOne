
import numpy as np
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import os, fire
import matplotlib.pyplot as plt

def chck_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def eval(latentX, cluster_num, sample_id, save_path, save_path_header):
    silhouette = []
    for cluster_num_ in cluster_num:
        predict_subtype = KMeans(n_clusters=cluster_num_, random_state=0).fit_predict(latentX)
        silhouette_ = silhouette_score(latentX, predict_subtype)
        silhouette.append(silhouette_)
        
        cluster_df = pd.DataFrame({"Sample": sample_id, "Subtype":predict_subtype} )
        save_cluster_path = save_path_header + f'_clusters={cluster_num_ }_clustering.csv'
        cluster_df.to_csv(save_cluster_path, header=True, index=False)
        
    fig, ax = plt.subplots()
    ax.plot(cluster_num, silhouette, marker="^", markersize=3, linewidth=1.5)
    ax.set_xlabel("#Cluster")
    ax.set_ylabel("Sihouette Score")
    plt.savefig(save_path)

    # cluster_num_wanted = cluster_num[ silhouette.index(max(silhouette ))  ]
    # predict_subtype = KMeans(n_clusters= cluster_num_wanted, random_state=0).fit_predict(latentX)
    # cluster_df = pd.DataFrame({"Sample": sample_id, "Subtype":predict_subtype} )
    # save_cluster_path = save_path_header + f'_clusters={cluster_num_wanted }_clustering.csv'
    # cluster_df.to_csv(save_cluster_path, header=True, index=False)

    return silhouette


def myeval(cluster_range="3-8"):
    alpha = [0.01, 0.02, 0.1, 1, 5]
    gamma = [0, 0.1, 1, 10, 100]
    low_dim = [2, 3, 4, 5, 6, 7]
    #low_dim = [7]
    c_min, c_max = cluster_range.split("-")
    cluster_num = [i for i in range(int(c_min), int(c_max) + 1 ) ]

    tdir = "./clusters/eval_cluster_num/"
    chck_dir(tdir)

    silhouette_summary = dict()
    silhouette_sm2 = pd.DataFrame(columns=["NMF_param","n_clusters","mean_siloutte_width"])
    for low_dim_ in low_dim:
        for alpha_ in alpha:
            for gamma_ in gamma:
                key = "lowDim=%d_alpha=%.2f_gamma=%.2f" % (low_dim_, alpha_, gamma_)

                latent_fi = "./NMF/lowDim=%d/" \
                            "X_lowDim=%d_alpha=%.2f_gamma=%.2f.csv" \
                            % (low_dim_, low_dim_, alpha_, gamma_)
                latentX_df = pd.read_csv(latent_fi, header=None, index_col=0)
                latentX = latentX_df.values
                sample_id = latentX_df.index

                save_path = tdir + "lowDim=%d_alpha=%.2f_gamma=%.2f_silhouette_score.png" \
                            % (low_dim_, alpha_, gamma_)
                save_path_header = tdir + "lowDim=%d_alpha=%.2f_gamma=%.2f" \
                            % (low_dim_, alpha_, gamma_)

                silhouette = eval(latentX, cluster_num, sample_id, save_path, save_path_header)

                silhouette_summary[key] = silhouette

                for sil, num in zip(silhouette, cluster_num):
                    sil_sm = pd.DataFrame({
                        "NMF_param": [ key ],
                        "n_clusters": [ num ],
                        "mean_siloutte_width": [ sil ]
                    })
                    silhouette_sm2 = silhouette_sm2.append(sil_sm, ignore_index=True)

    silhouette_df = pd.DataFrame.from_dict(silhouette_summary, orient="index",
                                        columns=["cluster_num="+str(num) for num in cluster_num] )
    silhouette_df_fout = tdir + "silhoutte_score_summary.csv"
    silhouette_df.to_csv(silhouette_df_fout, index=True, header=True)

    silhouette_fi_2 = tdir + "silhoutte_score_summary_2.csv"

    #df_sv = pd.read_csv("./record_log_rank_test_pvalue.csv")
    #silhouette_sm2 = silhouette_sm2.merge(df_sv, on="NMF_param", how="left")
    silhouette_sm2.to_csv(silhouette_fi_2, index=False)

if __name__ == '__main__':
    fire.Fire(myeval)