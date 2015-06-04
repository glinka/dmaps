#include <cmath>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
#include <string>
#include <sstream>
#include <util_fns.h>

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


/**
 * Creates swiss-roll dataset and returns the file location, relative to the dmaps directory
 * \param size size of the square to be sampled over
 * \param npts total number of points to sample from swissroll
 * \param r0 starting distance of swissroll's edge from origin
 * \note
 * Data is stored in csv file as:
 * Col1 | Col2 | Col3
 * -----------------
 *  x1  |  y1  |   z1
 *  x2  |  y2  |   z2
 */
std::string gen_swissroll(const int npts=1500, const int size=64, const double r0=0.25) {

  // if data already exists, do nothing
  std::stringstream ss;
  ss << "./inputdata/swissroll_" << npts << "_" << size << "_" << r0 << ".csv";
  std::string output_filename = ss.str();
  std::ifstream test_file_existence(output_filename);
  if(test_file_existence.good()) {
    test_file_existence.close();
    return output_filename;
  }
  else {
    test_file_existence.close();
  }
  // otherwise, generate new data and save
  // data is initially sampled over a square, then twisted into a three-dimensional swissroll shape
  std::mt19937 rng = gen_rng();
  std::vector<double> xvals(npts), yvals(npts), zvals(npts);
  // sample uniformly over grid
  for(int i = 0; i < npts; i++) {
    xvals[i] = size*gen_urn(rng);
    zvals[i] = size*gen_urn(rng);
  }
  // twist into swissroll
  for(int i = 0; i < npts; i++) {
    yvals[i] = sqrt(2*xvals[i] + r0*r0)*sin(sqrt(2*xvals[i] + r0*r0));
    xvals[i] = sqrt(2*xvals[i] + r0*r0)*cos(sqrt(2*xvals[i] + r0*r0));
  }

  // combine into one, save-able format
  std::vector< std::vector<double> > output_pts(npts, std::vector<double>(3));
  for(int i = 0; i < npts; i++) {
    for(int j = 0; j < 3; j++) {
      output_pts[i][0] = xvals[i];
      output_pts[i][1] = yvals[i];
      output_pts[i][2] = zvals[i];
    }
  }
  // save data
  util_fns::save_matrix(output_pts, output_filename);
  return output_filename;
}
