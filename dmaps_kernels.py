"""Stores alternative kernels to be used in the DMAPS algorithm

.. moduleauthor:: Alexander Holiday <holiday@alexanderholiday.com>

"""

import numpy as np

class objective_function_kernel:
    """A single-function class used to evaluate the modified DMAPS kernel between two points as motivated by Lafone's thesis. That is :math:`W_{ij}=exp(\\frac{\|x_i - x_j\|^2}{\epsilon} - \\frac{(of(x_i) - of(x_j))^2}{\epsilon^2})`

    Attributes:
        _epsilon (float): the DMAPS parameter :math:`\epsilon` to be used in kernel evaluations
    """

    def __init__(self, epsilon):
        # set epsilon
        self._epsilon = epsilon

    def __call__(self, pt1, pt2):
        """The function used to evaluate :math:`W_{ij}` between 'pt1' and 'pt2' with the prespecified value of :math:`\epsilon`. As a '__call__' method, the object itself is directly evaluated (see example below).

        .. note:
            pt[:-1] contains the parameter vector, while pt[-1] contains the objective function evaluation at that parameter set

        >>> of_kernel = objective_function_kernel(1e-1)
        >>> pt1 = np.arange(5); pt2 = np.arange(5,10)
        >>> print 'kernel evaluation between two pts:', of_kernel(pt1, pt2) # note how the object is called directly
        """
        return np.exp(-np.power(np.linalg.norm(pt1[:-1] - pt2[:-1]), 2)/self._epsilon - np.power(pt1[-1] - pt2[-1], 2)/np.power(self._epsilon, 2))

class gradient_kernel:
    """A single-function class used to evaluate the modified DMAPS kernel between two points as **defined** by Lafone's thesis. That is :math:`W_{ij}=exp(\\frac{\|x_i - x_j\|^2}{\epsilon} - \\frac{(<\nabla f_{x_i}, x_i-x_j>)^2}{\epsilon^2})`

    Attributes:
        _epsilon (float): the DMAPS parameter :math:`\epsilon` to be used in kernel evaluations
        _gradient (fn): gradient of :math:`f: \mathbb{R} \rightarrow \mathbb{R}^n`, accepting :math:`x \in \mathbb{R}^n` as argument


    """

    def __init__(self, epsilon, gradient):
        """Sets _epsilon and _gradient"""
        self._epsilon = epsilon
        self._gradient = gradient

    def __call__(self, pt1, pt2):
        """The function used to evaluate :math:`W_{ij}` between 'pt1' and 'pt2' with the prespecified value of :math:`\epsilon`. As a '__call__' method, the object itself is directly evaluated (see example below).

        >>> kernel = gradient_kernel(1e-1)
        >>> pt1 = np.arange(5); pt2 = np.arange(5,10)
        >>> print 'kernel evaluation between two pts:', kernel(pt1, pt2) # note how the object is called directly
        """
        return np.exp(-np.power(np.linalg.norm(pt1 - pt2), 2)/self._epsilon - np.power(np.dot(self._gradient(pt1), pt1 - pt2)/self._epsilon, 2))

class Data_Kernel:
    """Computes kernel between two points in parameter space, taking into account both the euclidean distance between parameters and the euclidean distance between model predictions at those parameters"""
    def __init__(self, epsilon, lam):
        self._epsilon = epsilon
        self._lam = lam

    def __call__(self, x1, x2):
        """Custom kernel given by: :math:`k(x_1, x_2) = e^{\frac{-1}{\lambda^2}(\frac{\| x_1 - x_2 \|^2}{\epsilon^2} + \|m(x_1) - m(x_2)\|^2)}` where :math:`m(x_i)` is the model prediction at parameter set :math:`x_i`

        Args:
            x1 (array): first data point in which x = [(parameters), (predictions)]
            x2 (array): second data point in which x = [(parameters), (predictions)]
        """
        return np.exp(-(np.power(np.linalg.norm(x1[0] - x2[0])/self._epsilon,2) + np.power(np.linalg.norm(x1[1] - x2[1]), 2))/(self._lam*self._lam))
