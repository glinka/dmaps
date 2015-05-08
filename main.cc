#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "dmaps_util_fns.h"
#include "dmaps.h"
#include "gaussian_kernel.h"
#include "gen_data.h"

int main(int argc, char *argv[]) {
  std::vector< std::string > possible_args;
  std::vector< std::string > input_filenames;
  possible_args.push_back("c");
  possible_args.push_back("create-folder");
  possible_args.push_back("s");
  possible_args.push_back("save");
  possible_args.push_back("m");
  // loop through all args and parse
  std::vector< std::string > output_data;
  std::string dir;
  int i = 1;
  int m = -1;
  // NON-OPT ARGS MUST PRECEED OPTIONAL ARGS
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
	else if(*arg == "s" || *arg == "save") {
	  // arg values will represent which data to be stored
	  // must be one or more of the following:
	  // "W", "eigvals", "eigvects"
	  // defaults to saving nothing
	  int j = 1;
	  while((i+j < argc) && (argv[i+j][0] != '-')) {
	    output_data.push_back(std::string(argv[i+j]));
	    j++;
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
    dir = "datadefault/";
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
    else if(*val == "eigs") {
      save_flags.SAVE_EIGVECTS = true;
      save_flags.SAVE_EIGVALS = true;
    }
    else {
      std::cout << "unrecognized output data entry" << std::endl;
      exit(1);
    }
  }
  // compute the mapping
  if(input_filenames.empty())
    input_filenames.push_back(gen_swissroll(120, 20));
  std::ifstream input_file(input_filenames[0]);
  std::vector< std::vector< double > > test_data = read_data(input_file);
  input_file.close();
  std::cout << "finished importing data from: " << input_filenames[0] << std::endl;
  int embedded_dim = 10;
  double weight_threshold = 1e-10;
  Gaussian_Kernel<double> gk_d(1e-4);
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // INIT THIS SHIT:
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  std::vector<double> eigvals;
  std::vector< std::vector<double> > eigvects, W;
  int dmaps_success = dmaps::map(test_data, gk_d, eigvals, eigvects, W, embedded_dim, weight_threshold);
  std::cout << "finished dmap " << std::endl;
  // create necessary files to store output data
  if(save_flags.SAVE_W) {
    std::stringstream ss;
    ss << dir << "W";
    std::ofstream W_output(ss.str());
    save_matrices(W_output, W);
    W_output.close();
    std::cout << "saved matrix W in: " << ss.str() << std::endl;
  }
  if(save_flags.SAVE_EIGVALS && save_flags.SAVE_EIGVECTS) {
    // save eigenvalues
    std::stringstream ss;
    ss << dir << "eigvals";
    std::ofstream eigvals_output(ss.str());
    save_vectors(eigvals_output, get_sorted_vals(eigvals), m);
    eigvals_output.close();
    std::cout << "saved eigenvalues in: " << ss.str() << std::endl;
    // save eigenvectors
    ss.str("");
    ss << dir << "eigvects";
    std::ofstream eigvects_output(ss.str());
    std::vector< int > sorted_indices = argsort(eigvals);
    save_matrices(eigvects_output, get_sorted_vectors(eigvects, sorted_indices), m);
    eigvects_output.close();
    std::cout << "saved eigenvectors in: " << ss.str() << std::endl;
  }
  else if(save_flags.SAVE_EIGVALS && !save_flags.SAVE_EIGVECTS) {
    std::stringstream ss;
    ss << dir << "eigvects";
    std::ofstream eigvects_output(ss.str());
    save_matrices(eigvects_output, eigvects);
    eigvects_output.close();
    std::cout << "saved eigenvectors in: " << ss.str() << std::endl;
  }
  else if(save_flags.SAVE_EIGVECTS && !save_flags.SAVE_EIGVALS) {
    std::stringstream ss;
    ss << dir << "eigvals";
    std::ofstream eigvals_output(ss.str());
    save_vectors(eigvals_output, eigvals);
    eigvals_output.close();
    std::cout << "saved eigenvalues in: " << ss.str() << std::endl;
  }
}
