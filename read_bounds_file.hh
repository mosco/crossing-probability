#ifndef __read_bounds_file_hh__
#define __read_bounds_file_hh__

#include <vector>
#include <string>
#include <utility>

std::pair<std::vector<double>, std::vector<double> > read_bounds_file(const std::string& filename);

#endif
