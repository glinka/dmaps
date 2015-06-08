#ifndef DMAPS_H
#define DMAPS_H
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <cmath>
#include "kernel_function.h"

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;

namespace dmaps {
  /**
   * Performs DMAP on 'input_data' using kernel 'kernel_fn', calculating a 'k'-dimensional embedding. <b> This method returns STL vectors and vector-of-vectors. </b>
   *
   * \tparam T input format of data, e.g. vector<double>, a matrix, some graph 
   * \param input_data vector in which each entry is a data point
   * \param kernel_fn kernel function used to evaluate each entry of 'W': \f$W_{ij} = k(x_i, x_j)\f$
   * \param eigvals preferably empty vector used to store DMAP's output eigenvalues
   * \param eigvects preferably empty vector of vectors used to store DMAP's output eigenvectors. Returned as size (k, n)
   * \note
   * Unlike most output from either Eigen or Python, <b> 'eigvects' contains the eigenvectors in its rows </b>
   * \param W the row-stochastic matrix used in the DMAP calculation
   * \param k dimension of embedding to compute, i.e. the number of eigval/eigvect pairs to find
   * \param weight_threshold threshold for elements of W, i.e. \f$W_{ij} < \f$ 'weight_threshold', set \f$W_{ij} = 0\f$
   * \returns indicator of success: 0 -> success, 1 -> failure
   */
  template <typename T>
    int map(const std::vector<T>& input_data, const Kernel_Function& kernel_fn, std::vector<double>& eigvals, std::vector< std::vector<double> >& eigvects, std::vector< std::vector<double> >& W, const int k=5, const double weight_threshold = 0);
    /**
   * Performs DMAP on 'input_data' using kernel 'kernel_fn', calculating a 'k'-dimensional embedding. <b> This method returns STL vectors and vector-of-vectors. This method does not return the W matrix. </b> As the entries of W are generally not of interest, this method simply saves copying data from Eigen::MatrixXd to an STL vector-of-vectors.
   *
   * \tparam T input format of data, e.g. vector<double>, a matrix, some graph 
   * \param input_data vector in which each entry is a data point
   * \param kernel_fn kernel function used to evaluate each entry of 'W': \f$W_{ij} = k(x_i, x_j)\f$
   * \param eigvals preferably empty vector used to store DMAP's output eigenvalues
   * \param eigvects preferably empty vector of vectors used to store DMAP's output eigenvectors. Returned as size (k, n)
   * \note
   * Unlike most output from either Eigen or Python, <b> 'eigvects' contains the eigenvectors in its rows </b>
   * \param k dimension of embedding to compute, i.e. the number of eigval/eigvect pairs to find
   * \param weight_threshold threshold for elements of W, i.e. \f$W_{ij} < \f$ 'weight_threshold', set \f$W_{ij} = 0\f$
   * \returns indicator of success: 0 -> success, 1 -> failure
   */
  template <typename T>
    int map(const std::vector<T>& input_data, const Kernel_Function& kernel_fn, std::vector<double>& eigvals, std::vector< std::vector<double> >& eigvects, const int k=5, const double weight_threshold = 0);
  /**
   * Performs DMAP on 'input_data' using kernel 'kernel_fn', calculating a 'k'-dimensional embedding. <b> This method returns Eigen-type vectors and matrices. </b>
   *
   * \tparam T input format of data, e.g. vector<double>, a matrix, some graph 
   * \param input_data vector in which each entry is a data point
   * \param kernel_fn kernel function used to evaluate each entry of 'W': \f$W_{ij} = k(x_i, x_j)\f$
   * \param eigvals preferably empty Eigen::VectorXd used to store DMAP's output eigenvalues
   * \param eigvects preferably empty Eigen::MatrixXd used to store DMAP's output eigenvectors. Returned as size (k, n)
   * \note
   * Unlike most output from either Eigen or Python, <b> 'eigvects' contains the eigenvectors in its rows </b>
   * \param W the row-stochastic matrix used in the DMAP calculation
   * \param k dimension of embedding to compute, i.e. the number of eigval/eigvect pairs to find
   * \param weight_threshold threshold for elements of W, i.e. \f$W_{ij} < \f$ 'weight_threshold', set \f$W_{ij} = 0\f$
   * \returns indicator of success: 0 -> success, 1 -> failure
   */
  template <typename T>
    int map(const std::vector<T>& input_data, const Kernel_Function& kernel_fn, Vector& eigvals, Matrix& eigvects, Matrix& W, const int k=5, const double weight_threshold = 0);
  /**
   * Investigates suitability of different DMAP kernels for use on 'input_data' by calculating \f$\sum_i \sum_j W_{ij}\f$ for each kernel. A kernel should be chosen from the region in which a plot of the output sum vs. kernel parameter is linear.
   *
   * \tparam T input format of data, e.g. vector<double>, a matrix, some graph 
   * \param input_data vector in which each entry is a data point
   * \param kernel_fns vector of kernel functions, each of which will be used to calculate a single value \f$\sum_i \sum_j k(x_i, x_j)\f$
   * \returns vector of \f$\sum_i \sum_j k(x_i, x_j)\f$ values for each input kernel
   */
  template <typename T>
    std::vector<double> test_kernels(const std::vector<T>& input_data, const std::vector<Kernel_Function>& kernel_fns);
}

#include "dmaps.tpp"

#endif
