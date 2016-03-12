import numpy as np
from joblib import Parallel, delayed

# Rows of X are patients/samples
# Cols of X are genes/features
def cormad(x, y):
    rad2 = np.sqrt(2)
    cost = 1.4826
    med_x = np.median(x)
    med_y = np.median(y)
    mad_x = cost * np.median(np.abs(x - med_x))
    mad_y = cost * np.median(np.abs(y - med_y))
    z_x = (x - med_x) / (cost * mad_x)
    z_y = (y - med_y) / (cost * mad_y)
    U = z_x + z_y
    V = z_x - z_y
    mad_U2 = ( cost * np.median( np.abs( U - np.median(U) )  ) )**2
    mad_V2 = ( cost * np.median( np.abs( V - np.median(V) )  ) )**2
    return (mad_U2 - mad_V2)  /  (mad_U2  + mad_V2)

def RobCor(X, n_jobs=1):
    n_samples, n_features = X.shape
    #corr = np.zeros((n_features, n_features))
    corr = []
    print "Estimating Correlations...",
    if n_jobs == 1:
        for i in range(n_features):
            for j in range(i):
                #corr[i, j] = cormad(X[:, i], X[:, j])
                corr.append(cormad(X[:, i], X[:, j]))
    else:
        corr = Parallel(n_jobs=n_jobs)(
                    delayed(cormad)(X[:, i], X[:, j]) for j in range(n_features) for i in range(n_features) if j < i
                    )
    print "Done."
    return corr
