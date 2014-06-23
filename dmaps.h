#ifndef DMAPS_H
#define DMAPS_H
#include <vector>

struct dmaps_output {
  std::vector< std::vector< double >* >* W;
  std::vector< std::vector< double >* >* eigvects;
  std::vector< double >* eigvals;
};
  
class dmaps {
 private:
 public:
  template <typename T>
    static dmaps_output* map(std::vector< T > &input_data, double (*kernel_fn)(const T&,const T&), const double weight_threshold = 0);
};

#include "dmaps.tpp"
#endif
