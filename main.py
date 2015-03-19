import numpy as np
import dmaps
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

def plot_xyz_data(x, y, z, color='b', **kwargs):
    """Plots three-dimensional data, used to display swissroll dataset"""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.grid(False)
    ax.scatter(x, y, z, c=color, **kwargs)
    ax.grid(False)
    plt.show()

def plot_plane(x, y, xlabel="", ylabel="", title="", color='b', **kwargs):
    """Plots two-dimensional data, used to display DMAPS results (two-dimensional embeddings)"""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ticksize = 24
    fontsize = 30
    plt.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.scatter(x, y, c=color, lw=0, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.show()

def gen_swissroll(n_thetas=60, n_zvals=10, var=0.1):
    """Generates a swissroll dataset in three dimensions

    Returns:
    swissroll (array): shape (n_zvals*n_thetas, 3) array in which each row represents the (x,y,z) coordinates of a point on the swissroll
    """
    theta_vals = np.linspace(0, 3*np.pi, n_thetas)
    r_vals = np.linspace(1, 3, n_thetas)
    x_vals = r_vals*np.cos(theta_vals)
    y_vals = r_vals*np.sin(theta_vals)
    z_vals = np.linspace(0, np.pi, n_zvals)
    swissroll = np.zeros([n_zvals*n_thetas, 3])
    for i in range(n_zvals):
        swissroll[i*n_thetas:(i+1)*n_thetas, 0] = x_vals + var*np.random.randn(n_thetas)
        swissroll[i*n_thetas:(i+1)*n_thetas, 1] = y_vals + var*np.random.randn(n_thetas)
        swissroll[i*n_thetas:(i+1)*n_thetas, 2] = z_vals[i] + var*np.random.randn(n_thetas)
    return swissroll

def dmaps_demo():
    """Demonstrates the DMAPS algorithm on a swissroll dataset using a predefined epsilon value"""
    data = gen_swissroll()
    epsilon = 1.25
    print 'Swissroll generated with', data.shape[0], 'points'
    print 'Displaying dataset'
    plot_xyz_data(data[:,0], data[:,1], data[:,2])
    start = time.clock()
    k = 4
    print 'Computing embedding'
    eigvals, eigvects = dmaps.embed_data(data, 5, epsilon=epsilon)
    print 'Lanczos solver took', str(time.clock() - start) + 's', 'to find top', k, 'eigenvectors'
    print 'Displaying dmaps embeddings'
    for i in range(1, k):
        for j in range(i+1, k):
            xlabel = r'$\Phi_' + str(i+1) + '$'
            ylabel = r'$\Phi_' + str(j+1) + '$'
            plot_plane(eigvals[i]*eigvects[:,i], eigvals[j]*eigvects[:,j], xlabel=xlabel, ylabel=ylabel, title='Embedding of swissroll with ' + xlabel + ' and ' + ylabel)

if __name__=="__main__":
    dmaps_demo()

