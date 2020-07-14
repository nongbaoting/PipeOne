import numpy as np
import pandas as pd
import os


def get_topk(X, weight, topk):
    idx = np.argsort(-weight)
    topk_idx = idx[:topk]
    return X.iloc[:, topk_idx]


proc_file_lst = ['alt_3prime_C3.confirmed.psi',
                 'alt_5prime_C3.confirmed.spi',
                 'apa', 'exon_skip_C3.confirmed.psi',
                 'fusion',
                 'intron_retention_C3.confirmed.psi',
                 'lncRNA',
                 'mult_exon_skip_C3.confirmed.psi',
                 'mutex_exons_C3.confirmed.psi',
                 'nonSilentGene',
                 'protein',
                 'retro',
                 'rnaEditing']

# best parameter settings
low_dim = 7
alpha = 5
gamma = 0.1

topk_ls = [50, 100, 200]

for topk in topk_ls:
    
    if not os.path.exists("./data/top%d/" % topk):
        os.makedirs("./data/top%d/" % topk)

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
