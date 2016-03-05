from cormad import cormad
import numpy as np
import time

M = np.loadtxt("data/RNASeq_samples_x_genes2.txt")

timeFirst = time.time()

cc = cormad(M)

elapsed = time.time()-timeFirst
print elapsed
print cc
