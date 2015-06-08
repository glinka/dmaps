#include "eigen_solvers.h"

namespace dmaps {


  template <typename T>
  int map(const std::vector<T>& input_data, const Kernel_Function& kernel_fn, Vector& eigvals, Matrix& eigvects, Matrix& W, const int k=5, const double weight_threshold = 0) {

    // calculate W entries
    int ndata = input_data.size();
    W = Matrix(ndata, ndata);
    for(int i = 0; i < ndata; i++) {
      for(int j = 0; j < ndata; j++) {
	W(i,j) = kernel_fn(input_data[i], input_data[j]);
      }
    }
    // ? for vectorization ?
    for(int i = 0; i < ndata; i++) {
      for(int j = 0; j < ndata; j++) {
	if(W(i,j) < weight_threshold) {
	  W(i,j) = 0;
	}
      }
    }

    // calculate row-stochastic matrix, calculate partial eigendecomp
    // normalize W to make symmetric
    Matrix D_half_inv = W.rowwise().sum().array().pow(-0.5).matrix().asDiagonal();
    Matrix S = D_half_inv*W*D_half_inv;
    Matrix V_ritz;
    Vector l_ritz;
    const int iram_maxiter=10*ndata, qr_maxiter=20*ndata;
    const int iram_success = eigen_solver::arnoldi_method_imprestart_hermitian(S, Vector::Ones(ndata), V_ritz, l_ritz, k, 2*k, iram_maxiter, qr_maxiter);
    eigvects = D_half_inv*V_ritz;
    eigvals = l_ritz;
    return iram_success == 1;
  }
    

  template <typename T>
  int map(const std::vector< T > &input_data, const Kernel_Function& kernel_fn, std::vector<double>& eigvals, std::vector< std::vector<double> >& eigvects, std::vector< std::vector<double> >& W, const int k, const double weight_threshold) {

    // set up "private" versions of output, to be converted into STL vectors later on
    Matrix _W, _eigvects;
    Vector _eigvals;
    int iram_success = map(input_data, kernel_fn, _eigvals, _eigvects, _W, k, weight_threshold);

    // copy data into output STL containers
    // unless _W is stored as row-major, I see now way of efficiently copying this data
    const int ndata = input_data.size();
    W = std::vector< std::vector<double> >(ndata, std::vector<double>(ndata));
    for(int i = 0; i < ndata; i++) {
      W[i][i] = _W(i,i);
      for(int j = i; j < ndata; j++) {
	W[i][j] = _W(i,j);
      }
    }
    double *eigvects_ptr = _eigvects.data();
    eigvects = std::vector< std::vector<double> >(k);
    for(int i = 0; i < k; i++) {
      eigvects[i] = std::vector<double>(eigvects_ptr + i*ndata, eigvects_ptr + (i+1)*ndata);
    }
    double *eigvals_ptr = _eigvals.data();
    eigvals = std::vector<double>(eigvals_ptr, eigvals_ptr + k);
    return iram_success == 1;
  }


  template <typename T>
  int map(const std::vector<T>& input_data, const Kernel_Function& kernel_fn, std::vector<double>& eigvals, std::vector< std::vector<double> >& eigvects, const int k=5, const double weight_threshold = 0) {

    // set up "private" versions of output, to be converted into STL vectors later on
    Matrix _W, _eigvects;
    Vector _eigvals;
    int iram_success = map(input_data, kernel_fn, _eigvals, _eigvects, _W, k, weight_threshold);

    // copy data into output STL containers
    double *eigvects_ptr = _eigvects.data();
    eigvects = std::vector< std::vector<double> >(k);
    const int ndata = input_data.size();
    for(int i = 0; i < k; i++) {
      eigvects[i] = std::vector<double>(eigvects_ptr + i*ndata, eigvects_ptr + (i+1)*ndata);
    }
    double *eigvals_ptr = _eigvals.data();
    eigvals = std::vector<double>(eigvals_ptr, eigvals_ptr + k);
    return iram_success == 1;
  }


  template <typename T>
  std::vector<double> test_kernels(const std::vector<T>& input_data, const std::vector<Kernel_Function>& kernel_fns) {
    const int npts = input_data.size();
    const int nkernels = kernel_fns.size();
    std::vector<double> w_sums(nkernels, 0);
    int sum_count = 0;
    // loop over each kernel
    for(auto kernel: kernel_fns) {
      // add up off-diagonal entries in upper right of matrix
      for(int i = 0; i < npts; i++) {
  	for(int j = i+1; j < npts; j++) {
  	  w_sums[sum_count] += kernel(input_data[i], input_data[j]);
  	}
      }
      // double to include off-diagonal in lower left
      w_sums[sum_count] *= 2;
      // finally, include diagonal elements
      for(int i = 0; i < npts; i++) {
	w_sums[sum_count] += kernel(input_data[i], input_data[i]);
      }
      sum_count++;
    }
    return w_sums;
  }


}
