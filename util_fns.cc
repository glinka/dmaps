#include "util_fns.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <fstream>

std::string create_directory(const std::string dir_basename) {
  int folder_counter = 0;
  std::string dir;
  bool isdir = true;
  std::stringstream ss;
  struct stat stat_dir;
  do {
    ss.str("");
    ss << dir_basename << folder_counter << "/";
    folder_counter++;
    dir = ss.str();
    int check = stat(dir.c_str(), &stat_dir);
    if(check == -1) {
      if(errno == ENOENT) {
	check = mkdir(dir.c_str(), 0700);
	if (check != 1) {
	  isdir = false;
	}
	else {
	  perror("mkdir error");
	  exit(1);
	}
      }
      else {
	perror("mkdir error");
	exit(1);
      }	  
    }
  } while (isdir);
  return dir;
}

std::vector< std::vector< double > > read_data(std::ifstream& input_file, const char delimiter) {
  //determine number of columns by reading first line
  std::string line;
  int column_count = 0;
  if(std::getline(input_file, line)) {
    std::stringstream ss(line);
    std::string temp;
    while(std::getline(ss, temp, delimiter)) {
      column_count++;
    }
  }
  else {
    std::cout << "empty input file" << std::endl;
    exit(1);
  }
  // move back to beginning of file
  input_file.seekg(0, std::ios::beg);
  std::string val;
  std::vector< std::vector< double > > output_data(column_count);
  int j;
  while(std::getline(input_file, line)) {
    std::stringstream ss(line);
    j = 0;
    while(std::getline(ss, val, delimiter)) {
      output_data[j++].push_back(atof(val.c_str()));
    }
  }
  // output_data[j++%current_column].push_back(val);
  if(j%column_count != 0) {
    std::cout << "missing data entries in input file" << std::endl;
    exit(1);
  }
  return output_data;
}

void save_matrices(std::ofstream& out_stream, const std::vector< std::vector< double >* >& data) {
  for(std::vector< std::vector< double >* >::const_iterator vect = data.begin(); vect != data.end(); vect++) {
    save_vectors(out_stream, *(*vect));
  }
}

void save_vectors(std::ofstream& out_stream, const std::vector< double >& data) {
  for(std::vector< double >::const_iterator val = data.begin(); val != (data.end()-1); val++) {
    out_stream << *val << " ";
  }
  out_stream << data.back() << std::endl;
}
