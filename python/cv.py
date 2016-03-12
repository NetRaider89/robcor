from cormad import RobCor
import numpy as np


def hardTresholding(X, ngrid=100, nsplits=100):
    R = RobCor(X)
    n_samples, n_features = X.shape
    tgrid = np.linspace(start=1e-5, stop=0.99, num=ngrid)
    n1 = n_samples - np.floor(n_samples / np.log(n_samples))
    FLOSSES = np.zeros((nsplits, ngrid))

    np.random.seed(19051977)

    for i in range(nsplits):
        idx = np.random.choice(a=range(n_samples), size=n1, replace=False)

        C1 = RobCor(X[idx, :])
        C2 = RobCor(np.delete(X, idx, axis=0))

        for k in range(ngrid):
            C1 = [0. if abs(c1) <= tgrid[k] else c1 for c1 in C1]
            FLOSSES[i, k] = np.linalg.norm(
                    [c1 - c2 for c1, c2 in zip(C1, C2)],
                    ord=2
                    )**2

    return (R, FLOSSES)

if __name__ == '__main__':
    data = np.loadtxt('data/RNASeq.txt').T

    R, FLOSSES = hardTresholding(data, 20, 20)
