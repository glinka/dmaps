import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

def plot_xy(x, y, xlabel="", ylabel="", title="", color='b', **kwargs):
    """Plots two-dimensional data, used to display DMAPS results (two-dimensional embeddings)"""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # hide ticks, too large, too busy
    plt.tick_params(axis='both', which='major', labelsize=0)
    ax.scatter(x, y, c=color, lw=0, **kwargs)
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
                plot_xy(np.power(eigvals[i], t)*eigvects[:,i], np.power(eigvals[j], t)*eigvects[:,j], xlabel=xlabel, ylabel=ylabel, s=50, **kwargs)
    if plot_3d:
        # loop through all the combinations
        for i in range(1, k):
            for j in range(i+1, k):
                for p in range(j+1, k):
                    xlabel = r'$\Phi_' + str(i+1) + '$'
                    ylabel = r'$\Phi_' + str(j+1) + '$'
                    zlabel = r'$\Phi_' + str(p+1) + '$'
                    plot_xyz(np.power(eigvals[i], t)*eigvects[:,i], np.power(eigvals[j], t)*eigvects[:,j], np.power(eigvals[p], t)*eigvects[:,p], xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, s=50, **kwargs)


