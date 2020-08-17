#ifndef __twosided_noncrossing_probability_hh__
#define __twosided_noncrossing_probability_hh__

#include <vector>

double ecdf_noncrossing_probability(int n, const std::vector<double>& b, const std::vector<double>& B, bool use_fft);

#endif
