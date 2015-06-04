#ifndef DMAPS_UTIL_FNS_H
#define DMAPS_UTIL_FNS_H
#include <vector>

namespace dmaps_utils {

  std::vector< int > argsort(const std::vector< double >& to_sort);
  std::vector< double > get_sorted_vals(const std::vector< double >& to_sort);
  std::vector< std::vector<double> > get_sorted_vectors(const std::vector< std::vector<double> >& to_sort, const std::vector< int >& sorted_indices);

}

#endif
