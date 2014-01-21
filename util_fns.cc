#include "util_fns.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <cstdio>

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
