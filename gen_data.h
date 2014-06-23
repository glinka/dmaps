#include <cmath>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
#include <string>
#include <sstream>

// pretty close to pi
const double PI = atan(1)*4;

std::mt19937 gen_rng() {
  // create a pseudo random number generator (Mersenne Twister algorithm)
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  return std::mt19937(seed);
}

double gen_urn(std::mt19937 &rng) {
  // get a uniform random number in the interval [0, 1)
  return rng()/((double) rng.max() + 1);
}


std::string gen_swissroll(const int nxy, const int nz, const double average_deviation = 0.35) {
  // creates swiss-roll dataset and
  // returns the file location,
  // relative to the dmaps directory
  /*
    Data is stored as
  
    x1 y1 z1
    x2 y2 z2
    ...
  
    as expected by read_data

  */
  // if data already exists, do nothing
  std::stringstream ss;
  ss << "inputdata/swissroll_" << nxy << "_" << nz << "_" << average_deviation << ".csv";
  std::string output_filename = ss.str();
  std::ifstream test_file_existence(output_filename);
  if(test_file_existence.good()) {
    test_file_existence.close();
    return output_filename;
  }
  else {
    test_file_existence.close();
  }
  std::mt19937 rng = gen_rng();
  const int n_xys = nxy, n_zs = nz;
  std::vector< double > xvals(n_xys), yvals(n_xys);
  double r = 1;
  const double r_final = 3;
  const double dr = (r_final - r) / n_xys;
  double theta = 0;
  const double theta_final = 3*PI;
  const double ds = (theta_final-theta)*(r_final)/n_xys;
  // calculate x and y values, in pairs
  // these pairs will then be broadcast across
  // all z values
  for(int i = 0; i < n_xys; i++) {
    xvals[i] = r*cos(theta);
    yvals[i] = r*sin(theta);
    r += dr;
    theta += ds/r;
  }
  double z = 0;
  const double z_final = 3;
  const double dz = (z_final - z) / n_zs;
  std::vector< double > zvals(n_zs);
  for(int i = 0; i < n_zs; i++) {
    zvals[i] = z;
    z += dz;
  }
  // create list of data vectors <x, y, z>
  // and add random noise to each entry
  std::vector< std::vector< double > > output_pts(n_xys*n_zs);
  int pt_counter = 0;
  for(int i = 0; i < n_xys; i++) {
    for(int j = 0; j < n_zs; j++) {
      output_pts[pt_counter].resize(3);
      output_pts[pt_counter][0] = xvals[i] + gen_urn(rng)*average_deviation;
      output_pts[pt_counter][1] = yvals[i] + gen_urn(rng)*average_deviation;
      output_pts[pt_counter][2] = zvals[j] + gen_urn(rng)*average_deviation;
      pt_counter++;
    }
  }
  // create file to save data
  std::ofstream swissroll(output_filename);
  for(std::vector< std::vector< double > >::const_iterator vec = output_pts.begin(); vec != output_pts.end(); vec++) {
    for(std::vector< double >::const_iterator val = (*vec).begin(); val != (*vec).end() - 1; val++) {
      swissroll << *val << ",";
    }
    swissroll << (*vec).back() << std::endl;
  }
  swissroll.close();
  return output_filename;
}
