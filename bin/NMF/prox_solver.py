import numpy as np
from numpy.linalg import norm


def prox_L1_solver(Y, alpha):
    nrow, ncol = Y.shape
    X = np.zeros((nrow, ncol))
    X[Y > alpha] = Y[Y > alpha] - alpha
    X[Y < -alpha] = Y[Y < -alpha] + alpha
    return X


def prox_L21_solver(Y, alpha, dim=0):
    nrow, ncol = Y.shape
    X = np.zeros((nrow, ncol))

    if dim == 0: # norm along rows
        for j in range(ncol):
            if norm(Y[:, j]) >= alpha:
                X[:, j] = (1 - alpha / norm(Y[:, j])) * Y[:, j]
    else:
        for i in range(nrow):
            if norm(Y[i, :]) >= alpha:
                X[i, :] = (1 - alpha/norm(Y[i, :])) * Y[i, :]

    return X


def L21norm(X, dim=0):
    # dim=0, norm along rows
    L21 = np.sum(np.apply_along_axis(norm, dim, X))
    return L21

