#ifndef GAUSSIAN_KERNEL_FN_H
#define GAUSSIAN_KERNEL_FN_H
#include <cmath>
#include <vector>
#include "kernel_function.h"

template <typename T>
class Gaussian_Kernel : public Kernel_Function< std::vector<T> > {
 public:
 Gaussian_Kernel(const double epsilon):epsilon_(epsilon) {}
  ~Gaussian_Kernel() {}
  virtual double kernel(const std::vector<T>& x1, const std::vector<T>& x2) const {
    int n = x1.size();
    double norm = 0;
    for(int i = 0; i < n; i++) {
      norm += std::pow(x1[i] - x2[i], 2);
    }
    return std::exp(-norm/epsilon_);
  };
 private:
    const double epsilon_;
};

#endif
