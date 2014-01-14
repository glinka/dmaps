import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def dmaps(data, epsilon):
    """
    m is number of data pts
    n is dimension of each pt
    """
    m = data.shape[1]
    n = data.shape[0]
    A = np.zeros([m, m])
    for i in range(m):
        for j in range(m):
            A[i,j] = np.exp(-np.power(np.linalg.norm(data[:,i] - data[:,j]), 2)/epsilon)
    D = np.identity(m)*np.sum(A, 1)
    W = np.dot(np.linalg.inv(D), A)
    eigvals, eigvects = np.linalg.eig(W)
    # eigvals = np.absolute(eigvals)
    sorted_indices = np.argsort(eigvals)
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[:,sorted_indices]
    plot_data(data[0,:], data[1,:], data[2,:], eigvects[:,-4], cmap='jet')
    plot_data(data[0,:], data[1,:], data[2,:], eigvects[:,-2], cmap='jet')
    plot_data(data[0,:], data[1,:], data[2,:], eigvects[:,-3], cmap='jet')
    plot_plane(eigvects[:,-2], eigvects[:,-3], eigvects[:,-2])
    plot_plane(eigvects[:,-2], eigvects[:,-4], eigvects[:,-2])

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
    ax.set_xlabel(r'$\Psi_{2}$')
    ax.set_ylabel(r'$\Psi_{3}$')
    ax.scatter(x, y, c=color, **kwargs)
    plt.show()

if __name__=="__main__":
    var = 0.05
    epsilon = 0.25
    n_thetas = 70
    theta_vals = np.linspace(0, 3*np.pi, n_thetas)
    r_vals = np.linspace(1, 3, n_thetas)
    x_vals = r_vals*np.cos(theta_vals)
    y_vals = r_vals*np.sin(theta_vals)
    n_zvals = 20
    z_vals = np.linspace(0, 2*np.pi, n_zvals)
    # z_rotation_angle = np.pi/4.
    # rotation_matrix = np.array([[np.cos(z_rotation_angle), -np.sin(z_rotation_angle), 0],
#                               [np.sin(z_rotation_angle), np.cos(z_rotation_angle), 0],
#                               [0, 0, 1]])
    data = np.zeros([3, n_zvals*n_thetas])
    for i in range(n_zvals):
        data[0,i*n_thetas:(i+1)*n_thetas] = x_vals + var*np.random.randn(n_thetas)
        data[1,i*n_thetas:(i+1)*n_thetas] = y_vals + var*np.random.randn(n_thetas)
        data[2,i*n_thetas:(i+1)*n_thetas] = z_vals[i] + var*np.random.randn(n_thetas)
    # data = np.dot(np.transpose(rotation_matrix), data)
    dmaps(data, epsilon)
