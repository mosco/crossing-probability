#ifndef __read_boundaries_file_hh__
#define __read_boundaries_file_hh__

#include <vector>
#include <string>
#include <utility>

std::pair<std::vector<double>, std::vector<double> > read_and_check_boundaries_file(std::string filename);

#endif
