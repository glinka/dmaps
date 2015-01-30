#include "util_fns.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>

std::vector< std::string > parse_arg(const char* char_arg, const std::vector< std::string > possible_args) {
  //determine which possible_args are single-lettered
  std::vector< std::string > shortcut_args;
  std::vector< std::string > full_args;
  for(std::vector< std::string >::const_iterator s = possible_args.begin(); s != possible_args.end(); s++) {
    if((*s).length() > 1) {
      full_args.push_back(*s);
    }
    else {
      shortcut_args.push_back(*s);
    }
  }
  //if one of full_args is found, return it
  //else search for multiple of shortcut_args
  std::string str_arg(char_arg);
  std::vector< std::string > found_args;
  for(std::vector< std::string >::const_iterator s =  full_args.begin(); s != full_args.end(); s++) {
    if(str_arg == *s) {
      found_args.push_back(str_arg);
      return found_args;
    }
  }
  for(std::vector< std::string >::const_iterator s =  shortcut_args.begin(); s != shortcut_args.end(); s++) {
    if(str_arg.find(*s) != std::string::npos) {
      found_args.push_back(*s);
    }
  }
  return found_args;
}

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
  /*

    Saves each row of data as a column in output_data
    
    The data is expected to be of the form:
    x1 y1 z1 ...
    x2 y2 z2 ...
      ...
    thus each "column" of output_data contains all
    input data of one variable

  */
  // determine number of columns by reading first line
  // count rows, not columns
  std::string line;
  int nrows = 0;
  while(std::getline(input_file, line)) {
    nrows++;
  }
  // int column_count = 0;
  // if(std::getline(input_file, line)) {
  //   std::stringstream ss(line);
  //   std::string temp;
  //   while(std::getline(ss, temp, delimiter)) {
  //     column_count++;
  //   }
  // }
  // else {
  //   std::cout << "empty input file" << std::endl;
  //   exit(1);
  // }
  // move back to beginning of file
  // need to clear EOF flag first
  input_file.clear();
  input_file.seekg(0);
  std::string val;
  std::vector< std::vector< double > > output_data(nrows);
  int j = 0;
  while(std::getline(input_file, line)) {
    std::stringstream ss(line);
    while(std::getline(ss, val, delimiter)) {
      output_data[j].push_back(atof(val.c_str()));
    }
    j++;
  }
  // output_data[j++%current_column].push_back(val);
  // if(j%column_count != 0) {
  //   std::cout << "missing data entries in input file" << std::endl;
  //   exit(1);
  // }
  return output_data;
}

void save_matrices(std::ofstream& out_stream, const std::vector< std::vector< double >* >& data, int m) {
  if(m < 0) {
    m = data.size();
  }
  for(int i = 0; i < m; i++) {
    save_vectors(out_stream, *data[i], m);
  }
}

void save_matrices(std::ofstream& out_stream, const std::vector< std::vector< double > >& data, int m) {
  if(m < 0) {
    m = data.size();
  }
  for(int i = 0; i < m; i++) {
    // want m vectors, but all data.size()
    // components of those vectors, so
    // pass save_vectors "m = -1"
    save_vectors(out_stream, data[i], -1);
  }
}

void save_vectors(std::ofstream& out_stream, const std::vector< double >& data, int m) {
  if(m < 0) {
    m = data.size();
  }
  for(int i = 0; i < m-1; i++) {
    out_stream << data[i] << ",";
  }
  out_stream << data[m-1] << std::endl;
}

bool argsort_comp(const std::vector< double >& p1, const std::vector< double >& p2) {
  /*
    typically would be p1[1] < p2[1]
    but we wish to sort in
    descending order
  */
  return (p1[1] > p2[1]);
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
  // assume vectors are stored in rows in "vectos" input
  std::vector< std::vector< double > > sorted_vectors(vectors.size());
  int i = 0;
  for(std::vector< int >::const_iterator index = sorted_indices.begin(); index != sorted_indices.end(); index++) {
    sorted_vectors[i++] = vectors[*index];
  }
  return sorted_vectors;
}
  
