from cormad import RobCor
import numpy as np
import time

M = np.loadtxt("data/data.txt")

timeFirst = time.time()

cc = RobCor(M, n_jobs=2)

elapsed = time.time()-timeFirst
print elapsed
print cc
