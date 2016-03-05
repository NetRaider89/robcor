import numpy as np

# Rows of X are patients/samples
# Cols of X are genes/features
def cormad(X):
    rad2 = np.sqrt(2)
    cost = 1.4826 

    n_samples, n_features = X.shape
    med = np.zeros(n_features)
    mad = np.zeros(n_features)
    z = np.zeros((n_features, n_samples))
    corr = np.zeros((n_features, n_features))

    # Precompute 
    print "Precomputing medians...",
    for i in range(n_features):
        med[i] = np.median(X[:, i])
        mad[i] = cost * np.median(np.abs(X[:, i] - med[i]))
        z[i, :] = (X[:, i] - med[i]) / (cost * mad[i])
    print "Done."

    print "Estimating Correlations...",
    for i in range(n_features):
        for j in range(i):
            U = z[i, :] + z[j, :]
            V = z[i, :] - z[j, :]
            mad_U2 = ( cost * np.median( np.abs( U - np.median(U) )  ) )**2
            mad_V2 = ( cost * np.median( np.abs( V - np.median(V) )  ) )**2
            corr[i, j] = (mad_U2 - mad_V2)  /  (mad_U2  + mad_V2)
    print "Done."
    return corr
