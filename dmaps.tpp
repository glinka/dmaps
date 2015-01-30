#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include "eigen_solvers.h"


namespace dmaps {
  template <typename T>
  int map(const std::vector< T > &input_data, const Kernel_Function<T>& kernel_fn, std::vector<double>& eigvals, std::vector< std::vector<double> >& eigvects, std::vector< std::vector<double> >& W_out, const int k, const double weight_threshold) {
    int ndata = input_data.size();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> W(ndata, ndata);
    for(int i = 0; i < ndata; i++) {
      for(int j = 0; j < ndata; j++) {
	W(i,j) = kernel_fn.kernel(input_data[i], input_data[j]);
	if(W(i,j) < weight_threshold) {
	  W(i,j) = 0;
	}
      }
    }
    Eigen::VectorXd Degs = W.rowwise().sum();
    Eigen::MatrixXd D_half_inv = Eigen::MatrixXd::Zero(ndata, ndata);
    for(int i = 0; i < ndata; i++) {
      D_half_inv(i,i) = pow(Degs[i], -0.5);
    }
    // normalize W
    Eigen::MatrixXd S = D_half_inv*W*D_half_inv;
    Eigen::MatrixXd V_ritz;
    Eigen::VectorXd l_ritz;
    const int iram_maxiter=10*ndata, qr_maxiter=20*ndata;
    const int iram_success = eigen_solver::arnoldi_method_imprestart_hermitian(S, Eigen::VectorXd::Ones(ndata), V_ritz, l_ritz, k, 2*k, iram_maxiter, qr_maxiter);
    Eigen::MatrixXd evects = D_half_inv*V_ritz;
    double *eigvals_ptr = l_ritz.data();
    double *eigvects_ptr = evects.data();
    double *W_ptr = W.data();
    W_out = std::vector< std::vector<double> >(k);
    eigvects = std::vector< std::vector<double> >(k);
    for(int j = 0; j < k; j++) {
      W_out[j] = std::vector<double>(W_ptr + j*ndata, W_ptr + (j+1)*ndata);
      eigvects[j] = std::vector<double>(eigvects_ptr + j*ndata, eigvects_ptr + (j+1)*ndata);
    }
    eigvals = std::vector<double>(eigvals_ptr, eigvals_ptr + k);
    return iram_success == 1;
  }
}
