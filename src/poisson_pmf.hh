#ifndef __poisson_pmf__
#define __poisson_pmf__

#include <cmath>
#include <cassert>

// Computes the probability of a Poisson random variable with intensity lambda:
// Pr[Pois(lambda)=k] = e^-lambda * lambda^k / k!
inline double poisson_pmf(double lambda, int k)
{
    assert(k >= 0);
    assert(lambda >= 0.0);

    if (lambda == 0.0) {
        return k == 0 ? 1.0 : 0.0;
    }
    return std::exp(-lambda + k*std::log(lambda) - std::lgamma(k+1));
}


class PoissonPMFGenerator {
public:
    PoissonPMFGenerator(int max_k);
    ~PoissonPMFGenerator();
    double evaluate_pmf(double lambda, int k) const;
    // Fills pmf_array_ptr with the PMF of a Poisson random variable:
    //     Pr[Pois(lambda) = 0], ..., Pr[Pois(lambda) = k]
    // Returns the smallest integer N such that for all indices i >= N we have that due to double-precision rounding
    //     Pr[Pois(lambda) = i] = 0
    void compute_array(int k, double lambda); 
    const double* get_array() const {return pmf_array_ptr;}
private:
    int max_k;
    double* log_gamma_LUT;
    double* pmf_array_ptr;
};

#endif
