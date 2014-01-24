#include "dmaps.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

template <typename T>
dmaps_output* dmaps::map(std::vector< T > &input_data, double (*kernel_fn)(T, T)) {
  int ndata = input_data.size();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> W(ndata, ndata);
  for(int i = 0; i < ndata; i++) {
    for(int j = 0; j < ndata; j++) {
      W(i,j) = (*kernel_fn)(input_data[i], input_data[j]);
    }
  }
  Eigen::VectorXd Degs = W.rowwise().sum();
  Eigen::MatrixXd D_inv = Eigen::MatrixXd::Zero(ndata, ndata);
  for(int i = 0; i < ndata; i++) {
    D_inv(i,i) = 1.0/Degs[i];
  }
  // normalize W
  // row stochasticize W
  // markovinize W
  Eigen::MatrixXd A = D_inv*W;
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  Eigen::MatrixXd evals = es.eigenvalues().real();
  Eigen::MatrixXd evects = es.eigenvectors().real();
  double *eigvals_ptr = evals.data();
  double *eigvects_ptr = evects.data();
  double *W_ptr = W.data();
  // possible to create W, eigvals, eigvects on heap and reference directly
  // into vector? seems not
  std::vector< std::vector< double >* >* W_ = new std::vector< std::vector< double >* >(ndata);
  std::vector< std::vector< double >* >* eigvects = new std::vector< std::vector< double >* >(ndata);
  std::vector< double >* eigvals = new std::vector< double >(ndata);
  for(int j = 0; j < ndata; j++) {
    (*W_)[j] = new std::vector< double >(W_ptr, W_ptr + ndata);
    (*eigvects)[j] = new std::vector< double >(eigvects_ptr, eigvects_ptr + ndata);
    (*eigvals)[j] = *(eigvals_ptr+j);
  }
  // my_allocator myalloc = new my_allocator(eigvects_ptr);
  // std::vector< double >* eigvals = new std::vector< double >(myalloc);
  dmaps_output* out = new dmaps_output;
  (*out).W = W_;
  (*out).eigvects = eigvects;
  (*out).eigvals = eigvals;
  return out;
} 