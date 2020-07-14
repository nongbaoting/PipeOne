import snf
from network_enhancement import network_enhancement as ne
import numpy as np

def snf_cal(D):

    view_num = len(D)
    k_neighs = 20
    mu = 0.4
    iter = 20

    aff_mat_ls = snf.make_affinity(D, K=k_neighs, mu=mu)

    if view_num != len(aff_mat_ls):
        assert "the number of views check fails"

    # denoise for each network
    order = 2
    alpha = 0.7
    ne_aff_mat_ls = [ne(aff_mat_ls[i], order, k_neighs, alpha) for i in range(view_num)]

    # fuse networks
    fuesd_network = snf.snf(ne_aff_mat_ls[0:3], K=k_neighs, t=iter)
    ne_fused_network = ne(fuesd_network, order, k_neighs, alpha)

    D = np.diag(ne_fused_network.sum(axis=1))
    L = D - ne_fused_network

    return L



