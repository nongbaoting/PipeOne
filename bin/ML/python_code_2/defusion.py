import numpy as np
from prox_solver import *
from numpy.linalg import norm
from copy import deepcopy

class DeFusion:
    def __init__(self, low_dim, alpha, gamma, silence, set_seed=False):

        self.beta = 1
        self.lambda_ = 1
        self.nrep = 10
        if set_seed:
            self.seed = np.array([i * 10 for i in range(self.nrep)], dtype=int)
        else:
            self.seed = None
        self.step = 50
        #self.maxiter = 10 # test only, for time saving
        self.maxiter = 1000
        self.maxinner = 5
        self.silence = silence

        self.alpha = alpha
        self.gamma = gamma
        self.low_dim = low_dim

    def eval_obj(self, D, X, Z, E, L):
        view_num = len(D)
        alpha = self.alpha
        beta = self.beta
        gamma = self.gamma
        lambda_ = self.lambda_
        obj = gamma*np.trace(X.T @ L @ X) + lambda_*np.sum(np.abs(X))

        for i in range(view_num):
            obj = obj + norm(D[i] - X @ Z[i] - E[i], ord="fro") ** 2 \
                  + alpha * np.sum(np.abs(E[i])) + beta * L21norm(Z[i])

        return obj

    def solver(self, D, L):
        view_num = len(D)
        sample_num = D[0].shape[0]
        sample_dim = [D[i].shape[1] for i in range(view_num)]

        beta = self.beta
        lambda_ = self.lambda_
        maxiter = self.maxiter
        maxinner = self.maxinner
        # number of repetition
        nrep = self.nrep
        step = self.step
        silence = self.silence

        low_dim = self.low_dim
        alpha = self.alpha
        gamma = self.gamma

        # record
        convergence_obj = []
        convergence_run = np.zeros((nrep, ))

        obj_val = []

        seed = self.seed

        for rp in range(nrep):
            # random initialization
            # set random seed
            if seed is not None:
                np.random.seed(seed[rp])
            initX = np.random.uniform(0, 2, (sample_num, low_dim))
            initZ = []
            initE = []
            for k in range(view_num):
                # set random seed
                if seed is not None:
                    np.random.seed(seed[rp] + k + 1)
                initZ.append(np.random.uniform(0, 2, (low_dim, sample_dim[k])))
                initE.append(np.zeros((sample_num, sample_dim[k])))

            X = initX
            Z = deepcopy(initZ)
            E = deepcopy(initE)

            start_eval = 0
            old_eval = np.inf
            new_eval = 0
            run = 1
            objval_on_rep = []

            while (abs(old_eval - new_eval) > np.abs(start_eval - new_eval)*1e-2) \
                    and (run <= maxiter):
                # solver for X
                for xinner in range(maxinner):
                    gradX = 2 * gamma * L @ X
                    for i in range(view_num):
                        gradX = gradX - 2 * (D[i] - X @ Z[i] - E[i]) @ Z[i].T

                    ZZ = sum([Z[i] @ Z[i].T for i in range(view_num)])
                    deltaX = 1 / (2*(norm(ZZ, ord="fro") + norm(L, ord="fro")))
                    X = prox_L1_solver(X - deltaX * gradX, lambda_*deltaX)
                    X[X < 0] = 0

                # solver for Z
                for i in range(view_num):
                    deltaZi = 1 / (2*norm(X.T @ X, ord="fro"))
                    for zinner in range(maxinner):
                        gradZi = -2 * X.T @ (D[i] - X @ Z[i] - E[i])
                        Zi = prox_L21_solver(Z[i] - deltaZi*gradZi, beta*deltaZi)
                        Zi[Zi < 0] = 0
                        Z[i] = Zi

                # solver for E
                for i in range(view_num):
                    E[i] = prox_L1_solver(D[i] - X @ Z[i], alpha/2)

                if run == 1:
                    obj = self.eval_obj(D, X, Z, E, L)
                    new_eval = obj
                    start_eval = new_eval

                    if not silence:
                        print("run  \t|obj    \t|")
                        print("%d   \t|%f   \t" % (run, new_eval))

                if run != 1 and np.mod(run, step) == 0:
                    old_eval = new_eval
                    obj = self.eval_obj(D, X, Z, E, L)
                    new_eval = obj

                    if not silence:
                        print("%d   \t|%f   \t" % (run, new_eval))

                obj_on_rep = self.eval_obj(D, X, Z, E, L)
                objval_on_rep.append(obj_on_rep)
                run = run + 1

            obj_val.append(new_eval)
            convergence_obj.append(objval_on_rep)
            convergence_run[rp] = run - 1

            if not silence:
                print("%d   \t|%d   \t|%f   " % (rp, run, new_eval))

            if obj_val[-1] == min(obj_val):
                X_f = X
                Z_f = deepcopy(Z)
                E_f = deepcopy(E)

        convergence = {"obj": convergence_obj, "run": convergence_run}
        return X_f, Z_f, E_f, convergence









