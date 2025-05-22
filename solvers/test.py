import numpy as np
from solver_so import solver_toa_upgrade_720_sat
import time

iters = 1000
start = time.time()
for _ in range(iters):
    
    data = np.random.randn(32)  # or your actual data
    sols = solver_toa_upgrade_720_sat.solver_toa_upgrade_sat_720(data)

stop = time.time()

print(f'{(start - stop)/iters:.2e}')
print(sols.shape)   # (2, 12)
print(sols)         # complex result
