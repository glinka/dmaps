import numpy as np
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def get_data(filename, header_rows=1, **kwargs):
    path_to_file = os.path.realpath(filename)
    # f = open(path_to_file, "r")
    # params_str = f.readline()
    # params = get_header_data(params_str)
    # f.close()
    data = np.genfromtxt(path_to_file, delimiter=",", skip_header=header_rows, **kwargs)
    return data#, params

def plot_dmap():
    m = 3
    eigvals = get_data("datadefault/eigvals", header_rows=0)
    # eigenvectors are stored in rows, not columns!
    eigvects = get_data("datadefault/eigvects", header_rows=0)
    sorted_indices = np.argsort(eigvals)
    eigvals = np.sort(eigvals)
    eigvects = eigvects[sorted_indices, :]
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.scatter(eigvals[0]*eigvects[0,:], eigvals[1]*eigvects[1,:])
    ax1.scatter(eigvals[0]*eigvects[0,:], eigvals[2]*eigvects[2,:])
    plt.show()

if __name__=="__main__":
    import sys
    plot_input(sys.argv[1])
