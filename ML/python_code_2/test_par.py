import numpy as np
from joblib import Parallel, delayed


def mult(factor, val):
    print("factor x val = %.2f" % (factor * val))
    return factor * val

if __name__=="__main__":
    vals = [2, 3, 4, 5, 6]
    facotr = 2

    Parallel(n_jobs=2)(delayed(mult)(facotr, val) for val in vals)