import numpy as np
import pickle
import matplotlib.pyplot as plt
from sklearn.model_selection import ParameterGrid


def chech_convergence(convergence, fout, title_):
    convergence_obj = convergence["obj"]
    convergence_run = convergence["run"]

    if len(convergence_obj) != len(convergence_run):
        assert "dimension matching check fails!"

    nrep = len(convergence_obj)
    fig, ax = plt.subplots()

    for i in range(nrep):
        obj = convergence_obj[i]
        iter = np.arange(convergence_run[i]) + 1
        ax.plot(iter, obj)

    ax.set(xlabel="iter", ylabel="loss", title=title_)
    fig.savefig(fout)


# low_dim = np.arange(2, 4, 1)
low_dim = np.array([7])
alpha = np.array([0.01, 0.02]) #, 0.1, 1, 5])
gamma = np.array([0, 0.1]) #, 1， 10, 100])
params = {"low_dim": low_dim, "alpha": alpha, "gamma": gamma}
params_grid = list(ParameterGrid(params))

for params_ in params_grid:
    tdir = "./results/lowDim=%d/" % params_["low_dim"]
    fi = tdir + "lowDim=%d_alpha=%.2f_gamma=%.2f.pkl" \
        % (params_["low_dim"], params_["alpha"], params_["gamma"])

    fread = open(fi, "rb")
    res = pickle.load(fread)
    convergence = res["convergence"]
    title_ = "lowDim=%d_alpha=%.2f_gamma=%.2f" \
             % (params_["low_dim"], params_["alpha"], params_["gamma"])

    fout = tdir + "convergenece_lowDim=%d_alpha=%.2f_gamma=%.2f.png" \
        % (params_["low_dim"], params_["alpha"], params_["gamma"])

    chech_convergence(convergence, fout, title_)


