import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import time
import os

def get_data(filename, header_rows=0, delim=',', **kwargs):
    path_to_file = os.path.realpath(filename)
    params = []
    if header_rows > 0:
        f = open(path_to_file, "r")
        params_str = f.readline()
        params = get_header_data(params_str)
        f.close()
    data = np.genfromtxt(path_to_file, delimiter=delim, skip_header=header_rows, **kwargs)
    return data, params

def get_header_data(header_str):
    BEGIN = 0
    comma = 1
    #create dict from header, based on key=value format in csv
    params = {}
    while comma > 0:
        equals = header_str.find("=")
        comma = header_str.find(",")
        params[header_str[BEGIN:equals]] = float(header_str[equals+1:comma])
        header_str = header_str[comma+1:]
    params[header_str[BEGIN:equals]] = float(header_str[equals+1:comma])
    #make integer, may not work especially well
    for key in params:
        if(params[key] % 1 == 0):
            params[key] = int(params[key])
    return params

def dmaps(data, epsilon):
    """
    m is number of data pts
    n is dimension of each pt
    """
    m = data.shape[1]
    n = data.shape[0]
    W = np.zeros([m, m])
    for i in range(m):
        for j in range(m):
            W[i,j] = np.exp(-np.power(np.linalg.norm(data[:,i] - data[:,j]), 2)/epsilon)
    D = np.identity(m)*np.sum(W, 1)
    D_half_inv = np.zeros((m,m))
    for i in range(m):
        D_half_inv[i,i] = np.power(D[i,i], -0.5)
    S = np.dot(np.dot(D_half_inv, W), D_half_inv)
    eigvals, eigvects = np.linalg.eigh(S)
    eigvects = np.dot(D_half_inv, eigvects)
    sorted_indices = np.argsort(eigvals)
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[:,sorted_indices]
    # return [eigvals, eigvects]
    # with eigvects in columns
    # and both sorted such that
    # eigval1 > eigval2

    np.savetxt('datadefault/W_frompython.csv', W, delimiter=',')
    print 'saved W matrix in: datadefault/W_frompython.csv'

    np.savetxt('datadefault/eigvals_frompython.csv', np.sort(eigvals), delimiter=',')
    print 'saved sorted eigvals in: datadefault/eigvals_frompython.csv'

    return [eigvals, eigvects]

def dmaps_slow(data, epsilon):
    """
    m is number of data pts
    n is dimension of each pt
    """
    m = data.shape[1]
    n = data.shape[0]
    W = np.zeros([m, m])
    for i in range(m):
        for j in range(m):
            W[i,j] = np.exp(-np.power(np.linalg.norm(data[:,i] - data[:,j]), 2)/epsilon)
    D = np.identity(m)*np.sum(W, 1)
    D_inv = np.zeros((m,m))
    for i in range(m):
        D_inv[i,i] = 1.0/D[i,i]
    A = np.dot(D_inv, W)
    eigvals, eigvects = np.linalg.eig(A)
    sorted_indices = np.argsort(eigvals)
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[:,sorted_indices]
    # return [eigvals, eigvects]
    # with eigvects in columns
    # and both sorted such that
    # eigval1 > eigval2
    return [eigvals, eigvects]

def plot_data(x, y, z, color, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.grid(False)
    ax.scatter(x, y, z, c=color, **kwargs)
    ax.grid(False)
    plt.show()

def plot_plane(x, y, color, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ticksize = 24
    fontsize = 30
    plt.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.scatter(x, y, c=color, lw=0, **kwargs)
    plt.show()

def gen_swissroll(n_thetas=60, n_zvals=10, var=0.1):
    theta_vals = np.linspace(0, 3*np.pi, n_thetas)
    r_vals = np.linspace(1, 3, n_thetas)
    x_vals = r_vals*np.cos(theta_vals)
    y_vals = r_vals*np.sin(theta_vals)
    z_vals = np.linspace(0, np.pi, n_zvals)
    swissroll = np.zeros([3, n_zvals*n_thetas])
    for i in range(n_zvals):
        swissroll[0,i*n_thetas:(i+1)*n_thetas] = x_vals + var*np.random.randn(n_thetas)
        swissroll[1,i*n_thetas:(i+1)*n_thetas] = y_vals + var*np.random.randn(n_thetas)
        swissroll[2,i*n_thetas:(i+1)*n_thetas] = z_vals[i] + var*np.random.randn(n_thetas)
    np.savetxt('inputdata/frompython.csv', np.transpose(swissroll), delimiter=',')
    print 'saved input dmaps data in: inputdata/frompython.csv'
    return swissroll

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-data', '-indata', nargs=1, type=str, default=None)
    args = parser.parse_args()
    if(args.input_data is not None):
        data, params = get_data(args.input_data[0])
        data = np.transpose(data)
        print data.shape
    else:
        data = gen_swissroll()

    epsilon = 1.25
    print 'swissroll generated'
    # plot_data(data[0,:], data[1,:], data[2,:], 'c', alpha=0.7)
    start = time.clock()
    eigvals, eigvects = dmaps(data, epsilon)
    print 'self-adjoint solver took', str(time.clock() - start) + 's'
    maxindex = 5
    eigvects_to_plot = np.array([eigvects[:,-i] for i in range(1, maxindex)])
    nvects = maxindex - 1
    for i in range(nvects-1):
        for j in range(i+1, nvects):
            plot_plane(eigvects_to_plot[i,:], eigvects_to_plot[j,:], color=eigvects_to_plot[i,:])
            print 'i =', maxindex-1-i, 'j =', maxindex-1-j
    # plot_data(data[0,:], data[1,:], data[2,:], eigvects[:,-4], cmap='jet')
    # plot_data(data[0,:], data[1,:], data[2,:], eigvects[:,-2], cmap='jet')
    # plot_data(data[0,:], data[1,:], data[2,:], eigvects[:,-3], cmap='jet')

    # for time comparison, run the following
    # which uses a general eigen-solver

    # # start = time.clock()
    # # eigvalsS, eigvectsS = dmaps_slow(data, epsilon)
    # # print 'general solver took', str(time.clock() - start) + 's'
    # # eigvects_to_plotS = np.array([eigvectsS[:,-i] for i in range(1, maxindex)])
    # # for i in range(nvects-1):
    # #     for j in range(i+1, nvects):
    # #         plot_plane(eigvects_to_plotS[i,:], eigvects_to_plotS[j,:], color=eigvects_to_plotS[i,:])
