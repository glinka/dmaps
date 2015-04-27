import util_fns as uf
import numpy as np
import scipy.sparse.linalg as spla


def _l2_distance(vector1, vector2):
    return np.linalg.norm(vector1 - vector2)

def _compute_embedding(W, k):
    """Calculates a partial ('k'-dimensional) eigendecomposition of W by first transforming into a self-adjoint matrix and then using the Lanczos algorithm.

    Args:
    W (array): symmetric, shape (npts, npts) array in which W[i,j] is the DMAPS kernel evaluation for points i and j
    k (int): the number of eigenvectors and eigenvalues to compute

    Returns:
    Returns:
    eigvals (array): shape (k) vector with first 'k' eigenvectors of DMAPS embedding sorted from largest to smallest
    eigvects (array): shape ("number of data points", k) array with the k-dimensional DMAPS-embedding eigenvectors. eigvects[:,i] corresponds to the eigenvector of the :math:'i^{th}'-largest eigenvalue, eigval[i].
    """
    m = W.shape[0]
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
        ndists = m*(m-1)/2
        # calc average, divide by 2 because each distance is double counted in W
        # important to do by hand and not by boolean indexing as certain off-diagonal values of W may be zero to numerical precision
        epsilon = np.sum(W)/(2.0*ndists)
    elif epsilon is "median":
        epsilon = np.median(W[W > 0])
    W = np.exp(-np.power(W, 2)/epsilon)
    eigvals, eigvects = _compute_embedding(W, k)
    return eigvals, eigvects


def embed_data_customkernel(data, k, kernel):
    """Computes the 'k'-dimensional DMAPS embedding of 'data' using the function 'kernel' to evaluate the DMAPS kernel between points and 'epsilon' as the characteristic radius of the neighborhood of each point. **Typically 'embed_data' should be used which employs the default exponential kernel with a potentially customized metric between points.**

    Args:
    data (iterable): typically a shape ("number of data points", "dimension of data") array containing the data as row vectors, but could be a list in which the :math:'i^{th}' entry contains :math:'i^{th}' data point, e.g. an adjacency matrix
    kernel (function): the kernel to be used in conjunction with 'data', accepting calls like kernel(data[i], data[j])
    k (int): number of dimensions to embed into


    Returns:
    eigvals (array): shape (k) vector with first 'k' eigenvectors of DMAPS embedding sorted from largest to smallest
    eigvects (array): shape ("number of data points", k) array with the k-dimensional DMAPS-embedding eigenvectors. eigvects[:,i] corresponds to the eigenvector of the :math:'i^{th}'-largest eigenvalue, eigval[i].
    """
    # m is number of data pts, len should work in all cases
    m = len(data)
    W = np.empty([m, m])
    # first populate W with metrics
    for i in range(m):
        W[i,i] = kernel(data[i], data[i])
        for j in range(i+1, m):
            W[i,j] = kernel(data[i], data[j])
            W[j,i] = W[i,j]
    eigvals, eigvects = _compute_embedding(W, k)
    return eigvals, eigvects
    
def epsilon_plot(epsilons, data, fraction_kept=1):
    """Displays a logarithmic plot of :math:`\sum_{i,j} W_{ij}(\epsilon)` versus :math:`\epsilon` over the range of epsilons provided as the first argument. Reasonable :math:`\epsilon` values will fall in the linear range of this figure. Also plots the mean and median of the squared distances for comparison.
    
    Args:
    epsilons (array): epsilon values at which to calculate :math:`\sum_{i,j} W_{ij}(\epsilon)`. Should span many orders of magnitude to ensure the diagram includes the asymptotes at :math:`\epsilon \rightarrow 0` and :math:`\epsilon \rightarrow \infty`
    data (array): size (n, p) array where 'n' is the number of data points and 'p' is the dimension of each point
    """
    import matplotlib.pyplot as plt
    data = np.copy(uf.thin_array(data, frac_to_keep=fraction_kept))
    n = data.shape[0]
    nepsilons = epsilons.shape[0]
    w_sums = np.empty((nepsilons))
    # loop over epsilons and calculate sum at each value
    for epsilon in enumerate(epsilons):
        w_sum = 0
        for i in range(n):
            for j in range(i+1, n):
                w_sum = w_sum + np.exp(-np.power(_l2_distance(data[i], data[j]), 2)/epsilon[1])
        # include diagonal (add n) and lower half (multiply by 2) of W
        w_sums[epsilon[0]] = w_sum*2 + n
    # calc mean and median, store values in array M
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            M[i,j] = _l2_distance(data[i], data[j])
    median = np.median(M[M > 0])
    mean = np.sum(M)/(n*(n-1)/2)
    # plot results
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(epsilons, w_sums)
    ax.axvline(x=mean, c='r', label=r'$\epsilon_{mean}$')
    ax.axvline(x=median, c='g', label=r'$\epsilon_{median}$')
    ax.set_xlabel(r'$\epsilon$')
    ax.set_ylabel(r'$\sum W_{ij}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc=2)
    plt.show(fig)
