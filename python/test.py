from cormad import cormad
import numpy as np
import time

M = np.loadtxt("data.txt")

timeFirst = time.time()

for i in range(0,20000):
	cc = cormad(M[:,0],M[:,1])

elapsed = time.time()-timeFirst
print elapsed
print cc
