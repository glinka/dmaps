import numpy as np
import dmaps
import plot_dmaps
import time

def gen_swissroll(n_thetas=20, n_zvals=20, var=0.5):
    """Generates a swissroll dataset in three dimensions

    Returns:
        swissroll (array): shape (n_zvals*n_thetas, 3) array in which each row represents the (x,y,z) coordinates of a point on the swissroll
    """
    # sample points over square, twist into roll
    size = 64
    r0 = 0.25
    npts = 1500
    xvals = size*np.random.uniform(size=npts)
    zvals = size*np.random.uniform(size=npts)
    yvals = np.sqrt(2*xvals + r0*r0)*np.sin(np.sqrt(2*xvals + r0*r0))
    xvals = np.sqrt(2*xvals + r0*r0)*np.cos(np.sqrt(2*xvals + r0*r0))
    swissroll = np.empty((npts, 3))
    swissroll[:,0] = xvals
    swissroll[:,1] = yvals
    swissroll[:,2] = zvals
    return swissroll

def dmaps_demo():
    """Demonstrates the DMAPS algorithm on a swissroll dataset using a predefined epsilon value"""

    data = gen_swissroll()
    epsilon = np.sqrt(5.0)
    print 'Swissroll generated with', data.shape[0], 'points'
    print 'Displaying dataset'
    plot_dmaps.plot_xyz(data[:,0], data[:,1], data[:,2], color=np.linalg.norm(data[:,:2], axis=1), s=80)
    # investigate proper epsilon
    print 'Investigating effect of epsilon on embedding (may take some time)'
    plot_dmaps.epsilon_plot(np.logspace(-3, 3, 10), data)
    start = time.clock()
    k = 30
    print 'Computing embedding'
    eigvals, eigvects = dmaps.embed_data(data, k, epsilon=epsilon)

    np.savetxt('./eigvects.csv', eigvects, delimiter=',')
    np.savetxt('./eigvals.csv', eigvals, delimiter=',')
    np.savetxt('./data.csv', data, delimiter=',')

    print 'Lanczos solver took', str(time.clock() - start) + 's', 'to find top', k, 'eigenvectors'
    print 'Displaying dmaps embeddings'
    for i in range(1, k):
        for j in range(i+1, k):
            xlabel = r'$\Phi_' + str(i+1) + '$'
            ylabel = r'$\Phi_' + str(j+1) + '$'
            plot_dmaps.plot_xy(eigvals[i]*eigvects[:,i], eigvals[j]*eigvects[:,j], xlabel=xlabel, ylabel=ylabel, title='Embedding dataset with ' + xlabel + ' and ' + ylabel, color=np.linalg.norm(data[:,:2], axis=1), s=50, scatter=True, hide_ticks=True)

if __name__=="__main__":
    # print 'no'
    dmaps_demo()

