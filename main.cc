#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <util_fns.h>
#include "dmaps_util_fns.h"
#include "dmaps.h"
#include "kernel_function.h"
#include "gen_data.h"

/**
 * Demonstrates the DMAPS algorithm on the canonical swissroll dataset, saving the output eigenvectors and eigenvalues in './outputdata/'
 */

int main(int argc, char *argv[]) {

  /* either no arguments are supplied and a swissroll is generated with 'gen_swissroll()', the default, while the DMAP uses epsilon = 2.5 */
  /* or */
  /* argc == 5, and the input is read as */
  /* npts (for swissroll) = argv[1] */
  /* size (for swissroll) = argv[2] */
  /* r0 (for swissroll) = argv[3] */
  /* epsilon (for DMAP) = argv[4] */
  // generate data based on presence/lack of input arguments
  int npts, size;
  double r0, epsilon;
  std::string input_filename;
  if(argc > 1) {
    npts = std::stoi(argv[1]);
    size = std::stoi(argv[2]);
    r0 = std::stof(argv[3]);
    epsilon = std::stof(argv[4]);
    input_filename = gen_swissroll(npts, size, r0);
  }
  else {
    input_filename = gen_swissroll();
    epsilon = 2.5;
  }

  // read data from file
  std::ifstream input_file(input_filename);
  std::vector< std::vector< double > > test_data = util_fns::read_data(input_file);
  std::cout << "finished importing data from: " << input_filename << std::endl;
  std::cout << "computing dmap... " << std::endl;

  // compute the mapping
  int embedded_dim = 10;
  double weight_threshold = 1e-10;
  // set up kernel, use standard e^(d(x, y)/ epsilon^2)
  Kernel_Function gaussian_kernel(epsilon); // standard gaussian kernel
  std::vector<double> eigvals;
  std::vector< std::vector<double> > eigvects, W;
  int dmaps_success = dmaps::map(test_data, gaussian_kernel, eigvals, eigvects, W, embedded_dim, weight_threshold);
  std::cout << "finished dmap " << std::endl;

  // save eigenvalues
  std::string output_dir = "./outputdata/";
  util_fns::save_vector(dmaps_utils::get_sorted_vals(eigvals), output_dir + "eigvals.csv");
  // save eigenvectors
  std::vector< int > sorted_indices = dmaps_utils::argsort(eigvals);
  util_fns::save_matrix(dmaps_utils::get_sorted_vectors(eigvects, sorted_indices), output_dir + "eigvects.csv");
  std::cout << "saved eigenvalues and eigenvectors in: " << output_dir << std::endl;

  // test test_kernels function over log-spaced epsilons
  const int nkernels = 20;
  const int lower_exp = -3, upper_exp = 3;
  Vector epsilons = Vector::LinSpaced(nkernels, lower_exp, upper_exp);
  for (int i = 0; i < nkernels; i++) {
    epsilons[i] = pow(10, epsilons[i]);
  }
  std::vector<Kernel_Function> kernels;
  for (int i = 0; i < nkernels; i++) {
    kernels.push_back(Kernel_Function(epsilons[i]));
  }
  std::vector<double> w_sums = dmaps::test_kernels(test_data, kernels);
  // save output
  util_fns::save_vector(w_sums, "./outputdata/w_sums.csv");
  util_fns::save_vector(std::vector<double>(epsilons.data(), epsilons.data()+nkernels), "./outputdata/epsilons.csv");

  return dmaps_success;
}
