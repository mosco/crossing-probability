#ifndef __read_boundaries_file_hh__
#define __read_boundaries_file_hh__

#include <vector>
#include <string>
#include <utility>

std::pair<std::vector<double>, std::vector<double> > read_boundaries_file(const std::string& filename);
void verify_one_sided_boundary_is_valid(const std::vector<double>& steps);
void verify_two_sided_boundaries_are_valid(const std::vector<double>& lower_bound_steps, const std::vector<double>& upper_bound_steps);

#endif
