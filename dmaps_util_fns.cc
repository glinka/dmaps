#include "dmaps_util_fns.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>

namespace dmaps_utils {

  bool argsort_comp(const std::vector< double >& p1, const std::vector< double >& p2) {
    /*
      typically would be p1[0] < p2[0]
      but we wish to sort in
      descending order
    */
    return (p1[0] > p2[0]);
  }

  std::vector< int > argsort(const std::vector< double >& to_sort) {
    /*
      return indices corresponding to a
      sorting of to_sort in decreasing order
    */
    int n = to_sort.size();
    std::vector< std::vector< double > > extended_vect(n);
    int i = 0;
    for(std::vector< std::vector< double > >::iterator vect = extended_vect.begin(); vect != extended_vect.end(); vect++) {
      vect->push_back(to_sort[i]);
      vect->push_back(i++);
    }
    std::sort(extended_vect.begin(), extended_vect.end(), argsort_comp);
    std::vector< int > sorted_indices(n);
    for(i = 0; i < n; i++) {
      sorted_indices[i] = extended_vect[i][1];
    }
    return sorted_indices;
  }

  bool comp(const double& p1, const double& p2) {
    return (p1 > p2);
  }

  std::vector< double > get_sorted_vals(const std::vector< double >& vals) {
    std::vector< double > sorted_vals = vals;
    std::sort(sorted_vals.begin(), sorted_vals.end(), comp);
    return sorted_vals;
  }

  std::vector< std::vector< double > > get_sorted_vectors(const std::vector< std::vector<double> >& vectors, const std::vector< int >& sorted_indices) {
    // assume vectors are stored in rows in "vectors" input
    std::vector< std::vector< double > > sorted_vectors(vectors.size());
    int i = 0;
    for(int index : sorted_indices) {
      sorted_vectors[i++] = vectors[index];
    }
    return sorted_vectors;
  }
  
}
