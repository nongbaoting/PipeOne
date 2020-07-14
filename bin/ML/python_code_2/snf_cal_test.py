from snf_cal import snf_cal
import snf
from snf import datasets
from prox_solver import prox_L1_solver

sim_data = datasets.load_simdata()
L = snf_cal(sim_data.data)
X = prox_L1_solver(L, 0.005)
print(" ")

