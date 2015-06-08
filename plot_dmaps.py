import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import util_fns as uf

def plot_xyz(x, y, z, xlabel="x", ylabel="y", zlabel="z", color='b', **kwargs):
    """Plots three-dimensional data, used to display swissroll dataset"""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.scatter(x, y, z, c=color, **kwargs)
    ax.grid(False)
    # hide labels, too squashed
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    plt.tick_params(axis='both', which='major', labelsize=0)
    plt.show()

def plot_xy(x, y, xlabel="", ylabel="", title="", color='b', xscale='linear', yscale='linear', scatter=False, hide_ticks=False, **kwargs):
    """Plots two-dimensional data, used to display DMAPS results (two-dimensional embeddings)"""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if hide_ticks:
        # hide ticks, too large, too busy
        plt.tick_params(axis='both', which='major', labelsize=0)
    # default to plot
    if scatter:
        ax.scatter(x, y, c=color, lw=0, **kwargs)
    else:
        ax.plot(x, y, c=color, lw=1, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.show()

def plot_embeddings(eigvects, eigvals, k='all', t=0, plot_2d=True, plot_3d=False, **kwargs):
    """Plots the "k Choose 2" different 2d embeddings based on the top 'k' eigenvectors from DMAPS

    Args:
        eigvects (array): columns contain DMAPS eigenvectors used to embed data
        eigvals (array): DMAPS eigenvalues
        .. note::
            It is assumed that the eigenvectors and eigenvalues are sorted in order of decreasing eigenvalue magnitude, i.e. eigvals[i] >= eigvals[i+1] for all i.
        k (int, 'all'): either an integer corresponding to the number of eigenvectors to consider or 'all' which considers all combinations
        t (float): the time parameter in DMAPS, typically zero in our work
    """
    if k is 'all':
        k = eigvals.shape[0]
    if plot_2d:
        # loop through all the combinations
        for i in range(1, k):
            for j in range(i+1, k):
                xlabel = r'$\Phi_' + str(i+1) + '$'
                ylabel = r'$\Phi_' + str(j+1) + '$'
                plot_xy(np.power(eigvals[i], t)*eigvects[:,i], np.power(eigvals[j], t)*eigvects[:,j], xlabel=xlabel, ylabel=ylabel, s=50, scatter=True, hide_ticks=True, **kwargs)
    if plot_3d:
        # loop through all the combinations
        for i in range(1, k):
            for j in range(i+1, k):
                for p in range(j+1, k):
                    xlabel = r'$\Phi_' + str(i+1) + '$'
                    ylabel = r'$\Phi_' + str(j+1) + '$'
                    zlabel = r'$\Phi_' + str(p+1) + '$'
                    plot_xyz(np.power(eigvals[i], t)*eigvects[:,i], np.power(eigvals[j], t)*eigvects[:,j], np.power(eigvals[p], t)*eigvects[:,p], xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, s=50, **kwargs)


def epsilon_plot(epsilons, data, fraction_kept=1):
    """Displays a logarithmic plot of :math:`\sum_{i,j} W_{ij}(\epsilon)` versus :math:`\epsilon` over the range of epsilons provided as the first argument. Reasonable :math:`\epsilon` values will fall in the linear range of this figure. Also plots the mean and median of the squared distances for comparison.
    
    Args:
        epsilons (array): epsilon values at which to calculate :math:`\sum_{i,j} W_{ij}(\epsilon)`. Should span many orders of magnitude to ensure the diagram includes the asymptotes at :math:`\epsilon \\rightarrow 0` and :math:`\epsilon \\rightarrow \infty`
        data (array): size (n, p) array where 'n' is the number of data points and 'p' is the dimension of each point
    """
    import matplotlib.pyplot as plt
    data = np.copy(uf.thin_array(data, frac_to_keep=fraction_kept))
    n = data.shape[0]
    nepsilons = epsilons.shape[0]
    w_sums = np.empty((nepsilons))
    # loop over epsilons and calculate sum at each value
    for k, epsilon in enumerate(epsilons):
        w_sum = 0
        for i in range(n):
            for j in range(i+1, n):
                w_sum = w_sum + np.exp(-np.power(np.linalg.norm(data[i] - data[j]), 2)/(epsilon*epsilon))
        # include diagonal (add n) and lower half (multiply by 2) of W
        w_sums[k] = w_sum*2 + n
    # calc mean and median, store values in array M
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            M[i,j] = np.linalg.norm(data[i] - data[j])
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


def kernel_plot(kernels, params, data, fraction_kept=1):
    """Displays a logarithmic plot of :math:`\sum_{i,j} W_{ij}(\epsilon)` versus :math:`\epsilon` over the range of epsilons provided as the first argument. Reasonable :math:`\epsilon` values will fall in the linear range of this figure. Also plots the mean and median of the squared distances for comparison.
    
    Args:
        kernels (list): kernel functions used to calculate :math:`W_{ij} = k(pt_i, pt_j)`. Typically there should be some :math:`\epsilon` parameter in the kernel function that varies over many orders of magnitude.
        params (array): vector of length 'nkernels' containing the different values of the parameter of interest used to create the different 'kernels'. Typically a vector of :math:`\epsilon` values.
        data (array): size (n, p) array where 'n' is the number of data points and 'p' is the dimension of each point
    """
    import matplotlib.pyplot as plt
    data = np.copy(uf.thin_array(data, frac_to_keep=fraction_kept))
    n = data.shape[0]
    nkernels = len(kernels)
    w_sums = np.empty((nkernels))
    # loop over epsilons and calculate sum at each value
    for k, kernel in enumerate(kernels):
        w_sum = 0
        for i in range(n):
            for j in range(i, n):
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
    plt.show(fig)

