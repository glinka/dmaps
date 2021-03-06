"""A flexible implementation of the DMAPS dimensionality reduction algorithm in Python.

.. moduleauthor:: Alexander Holiday <holiday@alexanderholiday.com>

"""

import numpy as np
import scipy.sparse.linalg as spla
from scipy.spatial.distance import pdist

def _l2_distance(vector1, vector2):
    """Returns the l2 norm of vector1 - vector2: :math:`\sqrt{\sum_i (x_i - y_i)^2}`"""
    return np.linalg.norm(vector1 - vector2)


def _compute_embedding_laplace_beltrami(W, k, symmetric=True):
    """Calculates a partial ('k'-dimensional) eigendecomposition of W by first transforming into a self-adjoint matrix and then using the Lanczos algorithm. **Unlike '_compute_embedding', this method normalizes W by an estimate of the local probability density at each point in order to remove the influence of nonuniform sampling from the embedding.** In this way, the eigenvalues and eigenvectors should actually approximate the eigenvalues and eigenvectors of the heat operator on the manifold.

    Args:
        W (array): symmetric, shape (npts, npts) array in which W[i,j] is the DMAPS kernel evaluation for points i and j
        k (int): the number of eigenvectors and eigenvalues to compute
        symmetric (bool): indicates whether the Markov matrix is symmetric or not. During standard useage with the default kernel, this will be true allowing for accelerated numerics. **However, if using custom_kernel(), this property may not hold.**

    Returns:
        eigvals (array): shape (k) vector with first 'k' eigenvectors of DMAPS embedding sorted from largest to smallest
        eigvects (array): shape ("number of data points", k) array with the k-dimensional DMAPS-embedding eigenvectors. eigvects[:,i] corresponds to the eigenvector of the :math:`i^{th}`-largest eigenvalue, eigval[i].
    """
    m = W.shape[0]
    # diagonal matrix D, inverse, sqrt
    local_density_estimate_inv = np.identity(m)/np.sum(W,1)
    W = np.dot(np.dot(local_density_estimate_inv, W), local_density_estimate_inv)
    # transform into self-adjoint matrix and find partial eigendecomp of this transformed matrix
    D_half_inv = np.identity(m)/np.sqrt(np.sum(W,1))
    eigvals, eigvects = None, None
    if symmetric:
        eigvals, eigvects = spla.eigsh(np.dot(np.dot(D_half_inv, W), D_half_inv), k=k) # eigsh (eigs hermitian)
    else:
        eigvals, eigvects = spla.eigs(np.dot(np.dot(D_half_inv, W), D_half_inv), k=k) # eigs (plain eigs)
    # transform eigenvectors to match W
    eigvects = np.dot(D_half_inv, eigvects)
    # sort eigvals and corresponding eigvects from largest to smallest magnitude  (reverse order)
    sorted_indices = np.argsort(np.abs(eigvals))[::-1]

    eigvals = eigvals[sorted_indices]
    # also scale eigenvectors to norm one
    eigvects = eigvects[:, sorted_indices]/np.linalg.norm(eigvects[:, sorted_indices], axis=0)
    return eigvals, eigvects    
    

def _compute_embedding(W, k, symmetric=True):
    """Calculates a partial ('k'-dimensional) eigendecomposition of W by first transforming into a self-adjoint matrix and then using the Lanczos algorithm.

    Args:
        W (array): symmetric, shape (npts, npts) array in which W[i,j] is the DMAPS kernel evaluation for points i and j
        k (int): the number of eigenvectors and eigenvalues to compute
        symmetric (bool): indicates whether the Markov matrix is symmetric or not. During standard useage with the default kernel, this will be true allowing for accelerated numerics. **However, if using custom_kernel(), this property may not hold.**

    Returns:
        eigvals (array): shape (k) vector with first 'k' eigenvectors of DMAPS embedding sorted from largest to smallest
        eigvects (array): shape ("number of data points", k) array with the k-dimensional DMAPS-embedding eigenvectors. eigvects[:,i] corresponds to the eigenvector of the :math:`i^{th}`-largest eigenvalue, eigval[i].
    """
    m = W.shape[0]
    # diagonal matrix D, inverse, sqrt
    D_half_inv = np.identity(m)/np.sqrt(np.sum(W,1))
    # transform into self-adjoint matrix and find partial eigendecomp of this transformed matrix
    eigvals, eigvects = None, None
    if symmetric:
        eigvals, eigvects = spla.eigsh(np.dot(np.dot(D_half_inv, W), D_half_inv), k=k) # eigsh (eigs hermitian)
    else:
        eigvals, eigvects = spla.eigs(np.dot(np.dot(D_half_inv, W), D_half_inv), k=k) # eigs (plain eigs)
    # transform eigenvectors to match W
    eigvects = np.dot(D_half_inv, eigvects)
    # sort eigvals and corresponding eigvects from largest to smallest magnitude  (reverse order)
    sorted_indices = np.argsort(np.abs(eigvals))[::-1]

    eigvals = eigvals[sorted_indices]
    # also scale eigenvectors to norm one
    eigvects = eigvects[:, sorted_indices]/np.linalg.norm(eigvects[:, sorted_indices], axis=0)
    return eigvals, eigvects


def embed_data(data, k, metric=_l2_distance, epsilon='mean', embedding_method=_compute_embedding):
    """Computes the 'k'-dimensional DMAPS embedding of 'data' using the function 'metric' to compute distances between points and 'epsilon' as the characteristic radius of the neighborhood of each point

    Args:
        data (iterable): typically a shape ("number of data points", "dimension of data") array containing the data as row vectors, but could be a list in which the :math:`i^{th}` entry contains :math:`i^{th}` data point, e.g. an adjacency matrix
            metric (function): the distance measure to be used in conjunction with 'data', accepting calls like metric(data[i], data[j])
        k (int): number of dimensions to embed into
        epsilon (string, float): one of either "median", "mean" or a float. If "median" or "mean", the "median" or "mean" of the distances between all points is used as the epsilon value. If a float is given, this value is used.

    Returns:
        eigvals (array): shape (k) vector with first 'k' eigenvectors of DMAPS embedding sorted from largest to smallest
        eigvects (array): shape ("number of data points", k) array with the k-dimensional DMAPS-embedding eigenvectors. eigvects[:,i] corresponds to the eigenvector of the :math:`i^{th}`-largest eigenvalue, eigval[i].

    >>> from test_dmaps import gen_swissroll
    >>> swissroll_data = gen_swissroll()
    >>> k = 15; epsilon = 2.5
    >>> eigvals, eigvects = dmaps.embed_data(data, k, epsilon)
    >>> from plot_dmaps import plot_embeddings
    >>> plot_embeddings(eigvects, eigvals, k=3)
    """
    # m is number of data pts, len should work in all cases
    m = len(data)
    W = np.empty([m, m])
    # first populate W with metrics
    for i in xrange(m):
        W[i,i] = 0
        for j in xrange(i+1, m):
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
    W = np.exp(-np.power(W, 2)/(epsilon*epsilon))
    eigvals, eigvects = embedding_method(W, k)
    return eigvals, eigvects


def embed_data_customkernel(data, k, kernel, symmetric=False):
    """Computes the 'k'-dimensional DMAPS embedding of 'data' using the function 'kernel' to evaluate the DMAPS kernel between points and 'epsilon' as the characteristic radius of the neighborhood of each point. **Typically 'embed_data' should be used which employs the default exponential kernel with a potentially customized metric between points.**

    Args:
        data (iterable): typically a shape ("number of data points", "dimension of data") array containing the data as row vectors, but could be a list in which the :math:`i^{th}` entry contains :math:`i^{th}` data point, e.g. an adjacency matrix
        kernel (function): the kernel to be used in conjunction with 'data', accepting calls like kernel(data[i], data[j])
        k (int): number of dimensions to embed into


    Returns:
        eigvals (array): shape (k) vector with first 'k' eigenvectors of DMAPS embedding sorted from largest to smallest
        eigvects (array): shape ("number of data points", k) array with the k-dimensional DMAPS-embedding eigenvectors. eigvects[:,i] corresponds to the eigenvector of the :math:`i^{th}`-largest eigenvalue, eigval[i].
    """
    # m is number of data pts, len should work in all cases
    m = len(data)
    W = np.empty([m, m])
    # first populate W with metrics
    for i in xrange(m):
        for j in xrange(m):
            W[i,j] = kernel(data[i], data[j])

    print 'finished constructing kernel matrix'

    eigvals, eigvects = _compute_embedding(W, k, symmetric=symmetric)
    return eigvals, eigvects

    
def epsilon_plot(data, filename=False):
    """Displays a logarithmic plot of :math:`L(\epsilon)` versus :math:`\epsilon` where :math:`L(\epsilon)` denotes the number of pairwise distances less than :math:`\epsilon`. Thus for small :math:`\epsilon`, :math:`L(\epsilon) \rightarrow n`, while large values of :math:`\epsilon` yield  :math:`L(\epsilon) \rightarrow n^2`. Reasonable :math:`\epsilon` values will fall in the linear range of this figure. Also plots the mean and median of the squared distances for comparison.
    
    Args:
        data (array): size (n, p) array where 'n' is the number of data points and 'p' is the dimension of each point
        filename (bool): the filename to save the figure as. If left to default value of False, figure is not saved
    >>> from test_dmaps import gen_swissroll
    >>> swissroll_data = gen_swissroll()
    >>> epsilon_plot(swissroll_data)
    """
    n = data.shape[0]
    # calculate pairwise distances and sort
    sorted_dists = np.sort(pdist(data, 'euclidean'))
    # find maximum and minimum distances
    logeps_min = np.floor(np.log10(sorted_dists[0])).astype(int)
    logeps_max = np.ceil(np.log10(sorted_dists[-1])).astype(int)
    # make epsilon values of interest
    neps = 20
    Ls = n*np.ones(neps)
    logepss = np.logspace(logeps_min - 1, logeps_max + 1, neps)
    current_index = 0
    for i, eps in enumerate(logepss):
        while current_index < sorted_dists.shape[0] and sorted_dists[current_index] < eps:
            current_index += 1
        Ls[i] += 2*current_index
    median = np.median(sorted_dists)
    mean = np.average(sorted_dists)
    # plot results
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(logepss, Ls)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.axvline(x=mean, c='r', label=r'$\epsilon_{mean}$')
    ax.axvline(x=median, c='g', label=r'$\epsilon_{median}$')
    ax.set_xlabel(r'$\epsilon$')
    ax.set_ylabel(r'$L(\epsilon)$')
    ax.legend(loc=2)
    if filename is not False:
        plt.savefig(filename)
    plt.show(fig)


def kernel_plot(kernels, params, data, filename=False):
    """Displays a logarithmic plot of :math:`\sum_{i,j} W_{ij}(\epsilon)` versus :math:`\epsilon` over the range of epsilons provided as the first argument. Reasonable :math:`\epsilon` values will fall in the linear range of this figure. Also plots the mean and median of the squared distances for comparison.
    
    Args:
        kernels (list): kernel functions used to calculate :math:`W_{ij} = k(pt_i, pt_j)`. Typically there should be some :math:`\epsilon` parameter in the kernel function that varies over many orders of magnitude.
        params (array): vector of length 'nkernels' containing the different values of the parameter of interest used to create the different 'kernels'. Typically a vector of :math:`\epsilon` values.
        data (array): size (n, p) array where 'n' is the number of data points and 'p' is the dimension of each point
        filename (bool): the filename to save the figure as. If left to default value of False, figure is not saved

    >>> from test_dmaps import gen_swissroll
    >>> swissroll_data = gen_swissroll()
    >>> epsilons = np.logspace(-3, 3 10)
    >>> kernels = [objective_function_kernel(eps) for eps in epsilons]
    >>> kernel_plot(kernels, epsilons, swissroll_data)
    """
    import matplotlib.pyplot as plt
    n = len(data) # data.shape[0]
    nkernels = len(kernels)
    w_sums = np.empty((nkernels))
    # loop over epsilons and calculate sum at each value
    for k, kernel in enumerate(kernels):
        w_sum = 0
        for i in xrange(n):
            for j in xrange(n):
                w_sum = w_sum + kernel(data[i], data[j])
        w_sums[k] = w_sum
    # plot results
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(params, w_sums)
    ax.set_xlabel(r'$\epsilon$')
    ax.set_ylabel(r'$\sum W_{ij}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc=2)
    if filename is not False:
        plt.savefig(filename)
    plt.show(fig)
