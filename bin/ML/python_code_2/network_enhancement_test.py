from scipy.io import loadmat
from network_enhancement import *
import pandas as pd

raw_data = loadmat("Raw_butterfly_network.mat")
w_butterfly0 = raw_data["W_butterfly0"]
labels = raw_data["labels"]

order = 2
alpha = 0.9
k = int(np.minimum(20, np.ceil(w_butterfly0.shape[0]/10)))
w_butterfly_ne = network_enhancement(w_butterfly0, order, k, alpha)
df = pd.DataFrame(w_butterfly_ne)
df.to_csv("W_butterfly_ne.csv", header=None, index=None)
