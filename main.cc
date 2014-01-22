#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "util_fns.h"
#include "dmaps.h"
#include "kernel_fn.h"

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
  

int main(int argc, char *argv[]) {
  std::vector< std::string > possible_args;
  std::vector< std::string > input_filenames;
  possible_args.push_back("c");
  possible_args.push_back("create-folder");
  possible_args.push_back("save");
  possible_args.push_back("savedata");
  possible_args.push_back("m");
  // loop through all args and parse
  std::vector< std::string > output_data;
  std::string dir;
  int i = 1;
  int m = 10;
  //OPTIONAL ARGS MUST SUCCEED NON-OPT ARGS
  while(i < argc) {
    std::vector< std::string > found_args;
    // optional args will begin with a hypen, and may contain
    // multiple arguments, i.e. "ls -lh"
    if(argv[i][0] == '-') {
      found_args = parse_arg(++argv[i], possible_args);
      // run through optional args
      for(std::vector< std::string >::const_iterator arg =  found_args.begin(); arg != found_args.end(); arg++) {
	if(*arg == "c" || *arg == "create-folder") {
	  //create new folder
	  dir = create_directory("data");
	  i++;
	}
	else if(*arg == "savedata" || *arg == "save") {
	  // arg values will represent which data to be stored
	  // must be one or more of the following:
	  // "W", "eigvals", "eigvects"
	  // defaults to saving eigenvalues and eigenvectors
	  int j = 1;
	  while((argv[i+j][0] != '-') && (i+j < argc)) {
	    output_data.push_back(std::string(argv[i+j++]));
	  } 
	  i += j;
	}
	else if(*arg == "m" || *arg == "embedded-dim") {
	  // only save top "m" eigenvalues and eigenvectors
	  m = atoi(argv[++i]);
	  i++;
	}
      }
    }
    // non-optional args consist solely of input
    // file names, where the input data will be read from
    else {
      input_filenames.push_back(std::string(argv[i]));
      i++;
    }
  }
  // set default directory name
  if(dir.empty())
    dir = "datadefault";
  // check which output_data values are legitimate
  // and open a corresponding file
  std::vector< std::string > data_options;
  std::vector< std::ofstream* > data_files;
  struct flags {
    bool SAVE_W = false, SAVE_EIGVECTS = false, SAVE_EIGVALS = false;
  } save_flags;
  for(std::vector< std::string >::const_iterator val = output_data.begin(); val != output_data.end(); val++) {
    if(*val == "W") {
      save_flags.SAVE_W = true;
    }
    else if(*val == "eigvects") {
      save_flags.SAVE_EIGVECTS = true;
    }
    else if(*val == "eigvals") {
      save_flags.SAVE_EIGVALS = true;
    }
    else {
      std::cout << "unrecognized output data entry" << std::endl;
      exit(1);
    }
  }
  // create necessary files to store output data
  std::ifstream file("test.csv");
  std::vector< std::vector< double > > test_data = read_data(file);
  dmaps_output* out = dmaps::map(test_data, gaussian_kernel);
  if(save_flags.SAVE_W) {
    std::stringstream ss(dir);
    ss << "W";
    std::ofstream W_output(ss.str());
    save_matrices(W_output, *((*out).W)));
    W_output.close();
  }
  if(save_flags.SAVE_EIGVECTS) {
    std::stringstream ss(dir);
    ss << "eigvects";
    std::ofstream eigvects_output(ss.str());
    save_matrices(eigvects_output, *out.eigvects);
    eigvects_output.close();
  }
  if(save_flags.SAVE_EIGVALS) {
    std::stringstream ss(dir);
    ss << "eigvals";
    std::ofstream eigvals_output(ss.str());
    save_vectors(eigvals_output, *out.eigvals);
    eigvals_output.close();
  }
  delete out.W;
  delete out.eigvects;
  delete out.eigvals;
}
