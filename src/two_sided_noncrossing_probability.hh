#ifndef __twosided_noncrossing_probability_hh__
#define __twosided_noncrossing_probability_hh__

#include <vector>

// Compute the probability that a homogeneous Poisson process in [0,1] of a given intensity will stay within
// lower and upper boundary functions.
//
// intensity
//     The intensity of the Poisson process. i.e. the expectation of the process at 1.
//
// lower_bound_steps
//     Times at which the lower boundary function crosses an integer for the first time.
//     e.g. lower_bound_steps[3] is the first time in which the lower bound is >= 3.
//
// upper_bound_steps
//     Times at which the upper boundary function crosses an integer for the last time.
//     e.g. upper_bound_steps[7] is the last time in which the upper bound is <= 7.
//
// use_fft
//     If true, use the fast O(n^2 log n) algorithm based on the Fast Fourier Transform.
//     If false, the asymptotically slower O(n^3) algorithm by Khmaladze&Shinjikashvili [2001] is used.
//     Note In cases where the lower and upper boundary are close or have a small number of steps, the O(n^3)
//     algorithm may be faster.
//
// endpoint
//     If endpoint is -1, this procedure will compute the noncrossing probability Pr[ g(t) <= p(t) <= h(t) for all 0<=t<=1]
//     where g(t) is specified by lower_bound_steps, h(t) by upper_bound_steps and p(t) is the poisson process.
//
//     If endpoint is >= 0, the following probability is computed instead:
//         Pr[ g(t) <= p(t) <= h(t) for all t and p(t) == endpoint]
//     i.e. The probability that the process will not cross the boundaries, and finish at the given endpoint.
// 
double poisson_process_noncrossing_probability(double intensity, const std::vector<double>& lower_bound_steps, const std::vector<double>& upper_bound_steps, bool use_fft, int endpoint);

// Compute the probability that an empirical CDF F^(t) with n points in [0,1] will stay within lower and upper boundary functions:
//     Pr[g(t) <= F^(t) <= h(t) for all 0<=t<=1]
//     where F^(t) is the empirical CDF of n uniform samples. i.e.
//     F^(t) = (number of X_i <= t)/n   where X_1,...,X_n ~ U[0,1].
//
// n
//     Number of samples from which the empirical CDF is constructed.
//
// lower_bound_steps
//     Times at which the lower boundary function crosses an integer for the first time.
//     e.g. lower_bound_steps[3] is the first time in which the lower bound is >= 3.
//
// upper_bound_steps
//     Times at which the upper boundary function crosses an integer for the last time.
//     e.g. upper_bound_steps[7] is the last time in which the upper bound is <= 7.
//
// use_fft
//     If true, use the fast O(n^2 log n) algorithm based on the Fast Fourier Transform.
//     If false, the asymptotically slower O(n^3) algorithm by Khmaladze&Shinjikashvili [2001] is used.
//     Note In cases where the lower and upper boundary are close or have a small number of steps, the O(n^3)
//     algorithm may be faster.
//
double ecdf_noncrossing_probability(int n, const std::vector<double>& lower_bound_steps, const std::vector<double>& upper_bound_steps, bool use_fft);

#endif
