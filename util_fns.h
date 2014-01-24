#ifndef UTIL_FNS_H
#define UTIL_FNS_H
#include <string>
#include <vector>
#include <fstream>

std::string create_directory(const std::string dir_basename);
std::vector< std::vector< double > > read_data(std::ifstream &input_file, const char delimiter=' ');
void save_matrices(std::ofstream& out_stream, const std::vector< std::vector< double >* >& data);
void save_vectors(std::ofstream& out_stream, const std::vector< double >& data);
#endif
