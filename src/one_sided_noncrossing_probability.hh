#ifndef __one_sided_noncrossing_probabiity_hh__
#define __one_sided_noncrossing_probabiity_hh__

#include <vector>

double binomial_process_lower_noncrossing_probability(int n, const std::vector<double>& lower_bound_steps);
double binomial_process_upper_noncrossing_probability(int n, const std::vector<double>& upper_bound_steps);

#endif
