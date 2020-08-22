import numpy as np
import pandas as pd
import os, fire, re
from defusion import DeFusion
from snf_cal import snf_cal
from sklearn.model_selection import ParameterGrid
import pickle
from joblib import Parallel, delayed

def chck_dir(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)


def save_result(X, Z, E, convergence, sample_id, var_names, fout):
    dict_ = {"X": X, "Z": Z, "E": E, "convergence": convergence,
             "sample_id": sample_id, "vars": var_names}
    f = open(fout, "wb")
    pickle.dump(dict_, f)
    f.close()


def save_feature_weight(Z, var_names, fout):
    weight = np.sum(abs(Z.T), axis=1)
    df_weight = pd.DataFrame(weight, index=var_names)
    df_weight.to_csv(fout, header=False, index=True)


def run_defusion(low_dim, alpha, gamma, view):

    data_dir = "./data/proc/"

    view_num = len(view)

    D = []
    var_names = []
    sample_id = None
    for v in range(view_num):
        # format of file names for each data type
        fi = data_dir + "proc_%s.csv" % view[v]
        df = pd.read_csv(fi, header=0, index_col=0)
        if v == 0:
            sample_id = list(df.index)
        D.append(1.0 * df.values)
        var_names.append(list(df.columns))

    L = snf_cal(D)

    tdir = "./NMF/lowDim=%d/" % low_dim
    chck_dir(tdir)
    print("running lowDim=%d_alpha=%.2f_gamma=%.2f" % (low_dim, alpha, gamma))
    fout = tdir + "lowDim=%d_alpha=%.2f_gamma=%.2f.pkl" % (low_dim, alpha, gamma)

    defusion_ = DeFusion(low_dim, alpha, gamma, silence=True)
    X, Z, E, convergence = defusion_.solver(D, L)
    save_result(X, Z, E, convergence, sample_id, var_names, fout)

    dfX = pd.DataFrame(X, index=sample_id)
    xfout = tdir + "X_lowDim=%d_alpha=%.2f_gamma=%.2f.csv" % (low_dim, alpha, gamma)
    dfX.to_csv(xfout, header=False, index=True)

    if len(Z) != len(var_names):
        assert "view number check fails!"

    for v in range(view_num):
        zfout = tdir + "weight_%s_lowDim=%d_alpha=%.2f_gamma=%.2f.csv" \
                % (view[v], low_dim, alpha, gamma)
        save_feature_weight(Z[v], var_names[v], zfout)

def MYRUN_DF(threads ):
    low_dim = np.arange(2, 8, 1)
    # low_dim = np.array([7])
    alpha = np.array([0.01, 0.02, 0.1, 1, 5])
    gamma = np.array([0, 0.1, 1, 10, 100])
    params = {"low_dim": low_dim, "alpha": alpha, "gamma": gamma}
    params_grid = list(ParameterGrid(params))

    view = []
    data_dir = "./data/proc/"
    for entry in os.scandir(data_dir):
        if re.search('.csv$', entry.name):
            the_view = re.sub('^proc_', '', entry.name)
            the_view = re.sub('\.csv$', '', the_view)
            view.append(the_view)
    print(view)

    Parallel(n_jobs=threads)(delayed(run_defusion)(params_["low_dim"],
                                                    params_["alpha"],
                                                    params_["gamma"], view) for params_ in params_grid)



if __name__=="__main__":
    fire.Fire(MYRUN_DF)

    ''' 
    view = ['alt_3prime_C3.confirmed.psi', 'alt_5prime_C3.confirmed.spi', 'apa',
            'exon_skip_C3.confirmed.psi', 'intron_retention_C3.confirmed.psi',
            'lncRNA', 'mult_exon_skip_C3.confirmed.psi', 'mutex_exons_C3.confirmed.psi',
            'protein', 'retro', 'rnaEditing', 'fusion', 'nonSilentGene']
    # 确定最优参数时注释
    
    low_dim = 7
    alpha = 5
    gamma = 0.1
    run_defusion(low_dim, alpha, gamma, view)
    

    ##这部分代码只是为了确定参数alpha, gamma  ### 怎么确定？
    low_dim = np.arange(2, 8, 1)
    low_dim = np.array([7])
    alpha = np.array([0.01, 0.02, 0.1, 1, 5])
    gamma = np.array([0, 0.1, 1, 10, 100])
    params = {"low_dim": low_dim, "alpha": alpha, "gamma": gamma}
    params_grid = list(ParameterGrid(params))

    core_num = 48
    Parallel(n_jobs=core_num)(delayed(run_defusion)(params_["low_dim"],
                                                    params_["alpha"],
                                                    params_["gamma"], view) for params_ in params_grid)
    '''
    
