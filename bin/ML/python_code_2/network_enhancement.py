"""
We transform the source code of Network Enhancement written in
MATLAB to python. The original MATLAB code is available at http://snap.stanford.edu/ne/

Reference
Wang, B., Pourshafeie, A., Zitnik, M. et al. Network enhancement as a
general method to denoise weighted biological networks. Nat Commun 9, 3108 (2018).
https://doi.org/10.1038/s41467-018-05469-x
"""

import numpy as np


def sub2ind(arr_shape, rows, cols):
    return rows*arr_shape[1] + cols


def ne_dn(w, type):
    eps = np.spacing(1)
    w = w*w.shape[0]
    w = np.double(w)
    D = np.sum(np.abs(w), axis=1) + eps

    if type == "ave":
        D = 1 / D
        D = np.diag(D)
        wn = D @ w
    elif type == "gph":
        D = 1 / np.sqrt(D)
        D = np.diag(D)
        wn = D @ w @ D
    return wn


def dominateset(aff_mat, nr_of_knn):
    B = np.argsort(-aff_mat, axis=1)
    A = -np.sort(-aff_mat, axis=1)
    res = A[:, 0:nr_of_knn]
    inds = np.tile(np.arange(aff_mat.shape[0]).reshape(-1, 1), (1, nr_of_knn))
    locs = B[:, 0:nr_of_knn]

    PNN_matrix1 = np.zeros(aff_mat.shape)
    idx = sub2ind(aff_mat.shape, inds.flatten(), locs.flatten())
    tmp = PNN_matrix1.flatten()
    tmp[idx] = res.flatten()
    PNN_matrix1 = tmp.reshape(aff_mat.shape)
    PNN_matrix = (PNN_matrix1 + PNN_matrix1.T) / 2
    return PNN_matrix


def transitionfield(w):
    eps = np.spacing(1)
    zero_index = np.nonzero(np.sum(w, axis=1)==0)
    w = w*w.shape[0]
    w = ne_dn(w, "ave")
    w_ = np.sqrt(np.sum(np.abs(w), axis=0)+eps)
    w = w / np.tile(w_, (w.shape[0], 1))
    w = w @ w.T
    wnew = w
    wnew[zero_index, ] = 0
    wnew[:, zero_index] = 0
    return wnew


def network_enhancement(W_in, order, k, alpha):
    sample_num = W_in.shape[0]
    W_in1 = W_in * (1 - np.eye(sample_num))
    # zero_index = np.where(np.sum(np.abs(W_in1), axis=0) > 0)
    zero_index = np.nonzero(np.sum(np.abs(W_in1), axis=0) > 0)
    W0 = W_in[np.ix_(zero_index[0], zero_index[0])]
    W = ne_dn(W0, "ave")
    W = (W + W.T) / 2

    DD = np.sum(np.abs(W0), axis=0)

    if len(np.unique(np.reshape(W.T, (-1, )))) == 2: # row first in python by default
        P = W
    else:
        P = (dominateset(np.abs(W), np.minimum(k, W.shape[0]-1)))*np.sign(W)

    P = P + np.eye(P.shape[0]) + np.diag(np.sum(np.abs(P.T), axis=0))

    P = transitionfield(P)
    eps = np.spacing(1)
    D, U = np.linalg.eig(P)
    sort_D_idx = np.argsort(D)
    D = np.sort(D)
    U = U[:, sort_D_idx]
    d = np.real(D - eps)
    d = (1 - alpha)*d / (1 - alpha*d**order)
    # D = np.diag(np.real(d))
    D = np.diag(np.real(d))
    W = U @ D @ U.T
    W = W * (1 - np.eye(W.shape[0])) / np.tile((1-np.diag(W)).reshape(-1, 1), (1, W.shape[0]))
    D = np.diag(DD)
    W = D @ W
    W[W<0] = 0
    W = (W + W.T) / 2
    W_out = np.zeros(W_in.shape)
    W_out[np.ix_(zero_index[0], zero_index[0])] = W

    return W_out



