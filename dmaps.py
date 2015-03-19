import numpy as np
import scipy.sparse.linalg as spla


def _l2_distance(vector1, vector2):
    return np.linalg.norm(vector1 - vector2)

def embed_data(data, k, metric=_l2_distance, epsilon='mean'):
    """Computes the 'k'-dimensional DMAPS embedding of 'data' using the function 'metric' to compute distances between points and 'epsilon' as the characteristic radius of the neighborhood of each point

    Args:
    data (iterable): typically a shape ("number of data points", "dimension of data") array containing the data as row vectors, but could be a list in which the :math:'i^{th}' entry contains :math:'i^{th}' data point, e.g. an adjacency matrix
    metric (function): the distance measure to be used in conjunction with 'data', accepting calls like metric(data[i], data[j])
    k (int): number of dimensions to embed into
    epsilon (string, float): one of either "median", "mean" or a float. If "median" or "mean", the "median" or "mean" of the distances between all points is used as the epsilon value. If a float is given, this value is used.

    Returns:
    eigvals (array): shape (k) vector with first 'k' eigenvectors of DMAPS embedding sorted from largest to smallest
    eigvects (array): shape ("number of data points", k) array with the k-dimensional DMAPS-embedding eigenvectors. eigvects[:,i] corresponds to the eigenvector of the :math:'i^{th}'-largest eigenvalue, eigval[i].
    """
    # m is number of data pts, len should work in all cases
    m = len(data)
    W = np.empty([m, m])
    # first populate W with metrics
    for i in range(m):
        W[i,i] = 0
        for j in range(i+1, m):
            W[i,j] = metric(data[i], data[j])
            W[j,i] = W[i,j]
    if epsilon is "mean":
        # number of distances, "m choose 2"
        ndists = m*(m-1)/2.0
        # calc average, divide by 2 because each distance is double counted in W
        epsilon = np.sum(W)/(2*ndists)
        # TESTING
        print 'testing diff of numpy average and my own. a complex calculation:', np.average(W[W > 1e-16]) - epsilon
    elif epsilon is "median":
        epsilon = np.median(W[W > 0])
    W = np.exp(-np.power(W, 2)/epsilon)
    # diagonal matrix D, inverse, sqrt
    D_half_inv = np.identity(m)/np.sqrt(np.sum(W,1))
    # transform into self-adjoint matrix and find partial eigendecomp of this transformed matrix
    eigvals, eigvects = spla.eigsh(np.dot(np.dot(D_half_inv, W), D_half_inv), k=k)
    # transform eigenvectors to match W
    eigvects = np.dot(D_half_inv, eigvects)
    # sort eigvals and corresponding eigvects from largest to smallest magnitude  (reverse order)
    sorted_indices = np.argsort(np.abs(eigvals))
    sorted_indices = sorted_indices[::-1]
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[:, sorted_indices]
    return eigvals, eigvects
