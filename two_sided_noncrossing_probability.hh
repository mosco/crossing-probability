#ifndef __twosided_noncrossing_probability_hh__
#define __twosided_noncrossing_probability_hh__

#include <vector>

double poisson_process_noncrossing_probability(double intensity, const std::vector<double>& lower_bound_steps, const std::vector<double>& upper_bound_steps, bool use_fft, int endpoint);
double binomial_process_noncrossing_probability(int n, const std::vector<double>& lower_bound_steps, const std::vector<double>& upper_bound_steps, bool use_fft);

#endif
