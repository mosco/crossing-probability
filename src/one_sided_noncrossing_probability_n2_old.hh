#ifndef __one_sided_noncrossing_probabiity_hh__
#define __one_sided_noncrossing_probabiity_hh__

#include <vector>

double ecdf_lower_noncrossing_probability_n2_old(int n, const std::vector<double>& lower_bound_steps);
double ecdf_upper_noncrossing_probability_n2_old(int n, const std::vector<double>& upper_bound_steps);

#endif
