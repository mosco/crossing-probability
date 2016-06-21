#include <cmath>
#include <cassert>
#include <algorithm>
#include "poisson_pmf.hh"
#include "aligned_mem.hh"

using namespace std;

PoissonPMFGenerator::PoissonPMFGenerator(int max_n)
{
    this->max_n = max_n;
    log_gamma_LUT = allocate_aligned_doubles(max_n+1);
    for (int i = 0; i < max_n+1; ++i) {
        log_gamma_LUT[i] = lgamma(i);
    }
}

PoissonPMFGenerator::~PoissonPMFGenerator()
{
    free_aligned_mem(log_gamma_LUT);
}

void PoissonPMFGenerator::compute_pmf(int n, double lambda, double* buffer)
{
    assert(n <= max_n);
    assert(n > 0);

    if (lambda == 0.0) {
        fill(&buffer[0], &buffer[n], 0);
        buffer[0] = 1.0;
    } else {
        double log_lambda = log(lambda);
        for (int i = 0; i < n; ++i) {
            buffer[i] = exp(-lambda + i*log_lambda - log_gamma_LUT[i+1]);
        }
    }
}
