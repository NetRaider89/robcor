import sys
import numpy as np
from cormad import cormad

selected_gene = int(sys.argv[1])
data = np.loadtxt('data/RNASeq.txt')
n_samples, n_features = data.shape

correlations = []
for i in range(n_samples):
    correlations.append(cormad(data[selected_gene, :], data[i, :]))

np.savetxt('test/python_correlations.txt', correlations)
