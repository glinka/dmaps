#ifndef UTIL_FNS_H
#define UTIL_FNS_H
#include <string>
#include <vector>
#include <fstream>

std::vector< std::string > parse_arg(const char* char_arg, const std::vector< std::string > possible_args);
std::string create_directory(const std::string dir_basename);
std::vector< std::vector< double > > read_data(std::ifstream &input_file, const char delimiter=',');
void save_matrices(std::ofstream& out_stream, const std::vector< std::vector< double >* >& data, int m = -1);
void save_matrices(std::ofstream& out_stream, const std::vector< std::vector< double > >& data, int m = -1);
void save_vectors(std::ofstream& out_stream, const std::vector< double >& data, int m = -1);
std::vector< int > argsort(const std::vector< double >& to_sort);
std::vector< double > get_sorted_vals(const std::vector< double >& to_sort);
std::vector< std::vector<double> > get_sorted_vectors(const std::vector< std::vector< double >* >& to_sort, const std::vector< int >& sorted_indices);
#endif
