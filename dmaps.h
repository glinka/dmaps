#ifndef DMAPS_H
#define DMAPS_H
#include <vector>
#include "kernel_function.h"

namespace dmaps {
  template <typename T>
    int map(const std::vector<T> &input_data, const Kernel_Function<T>& kernel_fn, std::vector<double>& eigvals, std::vector< std::vector<double> >& eigvects, std::vector< std::vector<double> >& W, const int k=5, const double weight_threshold = 0);
}

#include "dmaps.tpp"

#endif
