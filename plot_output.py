import numpy as np
import matplotlib.pyplot as plt
import os

def get_data(filename, header_rows=1, delim=' ', **kwargs):
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

def plot_xy(x, y, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y, **kwargs)
    plt.show()

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
    parser.add_argument('--eigenvalues', '-eigvals', nargs=1, type=str, default=None)
    parser.add_argument('--eigenvectors', '-eigvects', nargs=1, type=str, default=None)
    parser.add_argument('--input-data', '-indata', nargs=1, type=str, default=None)
    args = parser.parse_args()
    if (args.eigenvalues and args.eigenvectors) is not None:
        eigvals, params = get_data(args.eigenvalues[0], header_rows=0, delim=',')
        eigvects, params = get_data(args.eigenvectors[0], header_rows=0, delim=',')

        np.savetxt('datadefault/eigvals_fromC.csv', np.sort(eigvals), delimiter=',')
        print 'saved sorted eigvals in: datadefault/eigvals_fromC.csv'

        sorted_indices = np.argsort(eigvals)
        eigvals = np.sort(eigvals)
        eigvects = normalize_rows(eigvects)
        eigvects = eigvects[sorted_indices,:]
        print eigvals[-4:]
        # highly questionable rescaling
        # so everything isn't ~0
        # eigvals = eigvals/np.min(eigvals)
        # plot embedding
        t = 0
        maxindex = 6
        eigvects_to_plot = np.array([eigvects[-i,:] for i in range(1, maxindex)])
        eigvals_to_plot = np.array([eigvals[-i] for i in range(1, maxindex)])
        nvects = maxindex - 1
        for i in range(nvects-1):
            for j in range(i+1, nvects):
                plot_xy(np.power(eigvals_to_plot[i], t)*eigvects_to_plot[i,:], np.power(eigvals_to_plot[j], t)*eigvects_to_plot[j,:], c=eigvects_to_plot[i,:], lw=0, alpha=0.7)
                print 'i =', nvects-i, 'j =', nvects-j
    # used to visualize the test data
    if args.input_data is not None:
        plot_input(args.input_data[0])
