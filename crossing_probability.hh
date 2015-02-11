#ifndef __crossing_probability_hh__
#define __crossing_probability_hh__

#include <vector>

double crossing_probability(const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds);
double crossing_probability_fft(const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds);

#endif
