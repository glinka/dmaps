#ifndef KERNEL_FN_H
#define KERNEL_FN_H
#include <cmath>
#include <vector>
#include <iostream>

double gaussian_kernel(const std::vector< double >& x1, const std::vector< double >& x2) {
  double epsilon = 1.25;
  if(x1.size() != x2.size()) {
    std::cout << "array dimensions do not match" << std::endl;
    exit(1);
  }
  int n = x1.size();
  double norm = 0;
  for(int i = 0; i < n; i++) {
    norm += pow(x1[i] - x2[i], 2);
  }
  return exp(-norm/epsilon);
}

#endif
