#ifndef __noncrossing_probability_hh__
#define __noncrossing_probability_hh__

#include <vector>

double binomial_process_noncrossing_probability(const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds);
double binomial_process_noncrossing_probability_fft(const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds);

#endif
