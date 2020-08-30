#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>
#include "poisson_pmf.hh"
#include "aligned_mem.hh"

using namespace std;

PoissonPMFGenerator::PoissonPMFGenerator(int max_k)
{
    assert(max_k > 0);

    this->max_k = max_k;
    log_gamma_LUT = allocate_aligned_doubles(max_k+2);
    for (int i = 0; i < max_k+2; ++i) {
        log_gamma_LUT[i] = lgamma(i);
    }
   pmf_array_ptr = allocate_aligned_doubles(max_k+1);
    for (int i = 0; i < max_k+1; ++i) {
        pmf_array_ptr[i] = 0;
    }
}

PoissonPMFGenerator::~PoissonPMFGenerator()
{
    free_aligned_mem(pmf_array_ptr);
    free_aligned_mem(log_gamma_LUT);
}

double PoissonPMFGenerator::evaluate_pmf(double lambda, int k) const
{
    assert(lambda >= 0);
    assert(k >= 0);

    if (lambda == 0) {
        return (k == 0) ? 1.0 : 0.0;
    }

    return exp(-lambda + k*log(lambda) - log_gamma_LUT[k+1]);
}

void PoissonPMFGenerator::compute_array(int k, double lambda)
{
    assert(k >= 0);
    assert(k <= max_k);

    if (lambda < 0) {
            throw runtime_error("Expecting lambda>0 in PoissonPMFGenerator::compute_array()");
    }
    if (lambda == 0) {
        pmf_array_ptr[0] = 1;
        for (int i = 1; i < k+1; ++i) {
            pmf_array_ptr[i] = 0;
        }
        return;
    }

    double log_lambda = log(lambda);
    for (int i = 0; i < k+1; ++i) {
        pmf_array_ptr[i] = exp(-lambda + i*log_lambda - log_gamma_LUT[i+1]);
    }
}

