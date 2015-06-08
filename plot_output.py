import numpy as np
import matplotlib.pyplot as plt
import os
import util_fns as uf

import plot_dmaps

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

def normalize_rows(data):
    """"
    when data is a numpy ndarray
    this iterates over and 
    normalizes each row
    """
    return np.array([x/np.linalg.norm(x) for x in data])

def plot_input(filename):
    from mpl_toolkits.mplot3d import Axes3D
    data, params = get_data(filename, header_rows=0, delim=',')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(data[:, 0], data[:, 1], data[:, 2], lw=0)
    plt.show()

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--eigenvalues', '-eigvals', type=str, default=None)
    parser.add_argument('--eigenvectors', '-eigvects', type=str, default=None)
    parser.add_argument('--input-data', '-indata', type=str, default=None)
    parser.add_argument('--epsilons', type=str, default=None)
    parser.add_argument('--kernel-sums', type=str, default=None)
    args = parser.parse_args()
    if (args.eigenvalues and args.eigenvectors) is not None:
        # should be pre-sorted from low to high
        eigvals = uf.get_data(args.eigenvalues, header_rows=0, delim=',')
        eigvects = uf.get_data(args.eigenvectors, header_rows=0, delim=',')
        plot_dmaps.plot_embeddings(eigvects.T, eigvals)
    # used to visualize the test data
    if args.input_data is not None:
        plot_input(args.input_data)
    if (args.epsilons and args.kernel_sums) is not None:
        epsilons = uf.get_data(args.epsilons, header_rows=0)
        kernel_sums = uf.get_data(args.kernel_sums, header_rows=0)
        print epsilons, kernel_sums
        plot_dmaps.plot_xy(epsilons, kernel_sums, xlabel=r"$\epsilon$", ylabel="$\sum W_{ij}$", xscale='log', yscale='log')
        
