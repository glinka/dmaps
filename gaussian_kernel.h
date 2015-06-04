#ifndef GAUSSIAN_KERNEL_FN_H
#define GAUSSIAN_KERNEL_FN_H
#include <cmath>
#include <vector>
#include <utility>
#include "kernel_function.h"

template <typename T>
class Gaussian_Kernel : public Kernel_Function< std::vector<T> > {
public:
  // constructor
  Gaussian_Kernel(const double epsilon):_epsilon(epsilon) {}
  // copy constructor
  Gaussian_Kernel(const Gaussian_Kernel& gk):_epsilon(gk._epsilon)  {}
  // move constructor
  Gaussian_Kernel(Gaussian_Kernel&& gk): _epsilon(std::move(gk._epsilon)) {}
  /* // assignment operator */
  /* Gaussian_Kernel& operator=(const Gaussian_Kernel& gk) { */
  /*   _epsilon = gk._epsilon; */
  /*   return *this; */
  /* } */
  ~Gaussian_Kernel() {}
  virtual double kernel(const std::vector<T>& x1, const std::vector<T>& x2) const {
    int n = x1.size();
    double norm = 0;
    for(int i = 0; i < n; i++) {
      norm += std::pow(x1[i] - x2[i], 2);
    }
    return std::exp(-norm/(_epsilon*_epsilon));
  };
private:
  const double _epsilon;
};

#endif
