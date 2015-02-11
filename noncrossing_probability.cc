#include <algorithm>
#include <cmath>
#include <cassert>
#include "fftwconvolver.hh"

using namespace std;

enum BoundType {LOWER, UPPER, END};

struct Bound {
    double location;
    BoundType tag;
};

bool operator<(Bound b0, Bound b1)
{
    return (b0.location < b1.location);
}

vector<Bound> join_all_bounds(const vector<double>& lower_bounds, const vector<double>& upper_bounds)
{
    vector<Bound> bounds;

    assert(lower_bounds.size() == upper_bounds.size());
    int n = lower_bounds.size();

    Bound b;
    for (int i = 0; i < n; ++i) {
        b.location = lower_bounds[i];
        b.tag = LOWER;
        bounds.push_back(b);

        b.location = upper_bounds[i];
        b.tag = UPPER;
        bounds.push_back(b);
    }

    sort(bounds.begin(), bounds.end());

    return bounds;
}

// inline __float128 poisson_pmf_float128(__float128 lambda, int k)
// {
//     assert(k >= 0);
//     assert(lambda >= 0.0);
// 
//     if (lambda == 0.0) {
//         return k == 0 ? 1.0 : 0.0;
//     }
//     __float128 log_pmf = -lambda + k*logq(lambda) - lgammaq(k+1);
//     return expq(log_pmf);
// }
//
// Computes the probability of a Poisson random variable with intensity lambda:
// Pr[Pois(lambda)=k] = e^-lambda * lambda^k / k!
inline double poisson_pmf(double lambda, int k)
{
    assert(k >= 0);
    assert(lambda >= 0.0);

    if (lambda == 0.0) {
        return k == 0 ? 1.0 : 0.0;
    }
    double log_pmf = -lambda + k*log(lambda) - lgamma(k+1);
    return exp(log_pmf);
}

// An efficient implementation of the algorithm for the computation of non-crossing probability for a binomial process given in the paper:
//     Khmaladze and Shinjikashvili (2001) "Calculation of noncrossing probabilities for Poisson processes and its corollaries"
double binomial_process_noncrossing_probability(const vector<double>& lower_bounds, const vector<double>& upper_bounds)
{
    assert(lower_bounds.size() == upper_bounds.size());
    int n = lower_bounds.size();

    vector<Bound> bounds = join_all_bounds(lower_bounds, upper_bounds);
    Bound b;
    b.location = 1.0;
    b.tag = END;
    bounds.push_back(b);

    vector<double> Qs0(n+1, 0.0);
    Qs0[0] = 1.0;
    vector<double> Qs1(n+1, 0.0);

    vector<double>* buffers[] = {&Qs0, &Qs1};

    double prev_location = 0.0;
    int lower_bound_count = 0;
    int upper_bound_count = 0;

    vector<double> pmf(n+1, 0.0);

    for (unsigned int i = 0; i < bounds.size(); ++i) {
        // cout << "--------------------------------------------\n";
        // cout << "Iteration " << i << endl;

        const vector<double>& from = *buffers[i % 2];
        vector<double>& to = *buffers[(i+1) % 2];

        // cout << "From:\n";
        // print_vector(from);
        // cout << endl;

        // cout << "lower_bound_count: " << lower_bound_count << endl;
        // cout << "upper_bound_count: " << upper_bound_count << endl;

        double location = bounds[i].location;

        for (int j = 0; j < lower_bound_count - upper_bound_count + 1; ++j) {
            pmf[j] = poisson_pmf(n*(location-prev_location), j);
        }

        for (int j = upper_bound_count; j <= lower_bound_count; ++j) {
            double convolution_at_j = 0.0;
            for (int k = upper_bound_count; k <= j; ++k) {
                //convolution_at_j += from[k] * poisson_pmf(n*(location-prev_location), j-k);
                convolution_at_j += from[k] * pmf[j-k];
            }
            to[j] = convolution_at_j;
        }
        // cout << "convolution: ";
        // print_vector(to);

        //cout << i << ": ";
        //print_double_array(&to[upper_bound_count], lower_bound_count-upper_bound_count+1);

        BoundType tag = bounds[i].tag;
        if (tag == LOWER) {
            ++lower_bound_count;
        } else if (tag == UPPER) {
            ++upper_bound_count;
        } else {
            assert(tag == END);
        }
        prev_location = location;

        //cout << "To:\n";
        //print_vector(to);
        //cout << endl;
    }

    vector<double>& last_to_buffer = *buffers[bounds.size() % 2];
    double poisson_process_noncrossing_probability = last_to_buffer[n];
    return poisson_process_noncrossing_probability / poisson_pmf(n, n);
}


double binomial_process_noncrossing_probability_fft(const vector<double>& lower_bounds, const vector<double>& upper_bounds)
{
    assert(lower_bounds.size() == upper_bounds.size());
    int n = lower_bounds.size();

    // cout << "lower_bounds: ";
    // print_vector(lower_bounds);
    // cout << "upper_bounds: ";
    // print_vector(upper_bounds);

    vector<Bound> bounds = join_all_bounds(lower_bounds, upper_bounds);
    Bound b;
    b.location = 1.0;
    b.tag = END;
    bounds.push_back(b);

    vector<double> Qs0(n+1, 0.0);
    Qs0[0] = 1.0;
    vector<double> Qs1(n+1, 0.0);

    vector<double>* buffers[] = {&Qs0, &Qs1};

    double prev_location = 0.0;
    int lower_bound_count = 0;
    int upper_bound_count = 0;

    //vector<typename fftw_traits<T>::complex_type> tmp(2*n);
    FFTWConvolver convolver(n+1);
    for (unsigned int i = 0; i < bounds.size(); ++i) {
        // cout << "--------------------------------------------\n";
        // cout << "Iteration " << i << endl;

        const vector<double>& from = *buffers[i % 2];
        vector<double>& to = *buffers[(i+1) % 2];

        //cout << "From:\n";
        //print_array(&from[0], n+1);
        //cout << endl;

        //cout << "lower_bound_count: " << lower_bound_count << endl;
        //cout << "upper_bound_count: " << upper_bound_count << endl;

        double location = bounds[i].location;

        int cur_size = lower_bound_count - upper_bound_count + 1;
        vector<double> pmf(cur_size);
        for (unsigned int j = 0; j < pmf.size(); ++j) {
            pmf[j] = poisson_pmf(n*(location-prev_location), j);
        }
        //cout << "pmf: ";
        //print_array(&pmf[0], n+1);
        convolver.convolve_same_size(cur_size, &pmf[0], &from[upper_bound_count], &to[upper_bound_count]);
        //cout << "convolution result: ";
        //print_array(&to[0], n+1);

        //cout << i << ": ";
        //print_float128_array(&to[upper_bound_count], cur_size);

        BoundType tag = bounds[i].tag;
        if (tag == LOWER) {
            ++lower_bound_count;
        } else if (tag == UPPER) {
            ++upper_bound_count;
        } else {
            assert(tag == END);
        }
        prev_location = location;

        //cout << "To:\n";
        //print_vector(to);
        //cout << endl;
    }

    vector<double>& last_to_buffer = *buffers[bounds.size() % 2];
    double poisson_process_noncrossing_probability = last_to_buffer[n];

    return poisson_process_noncrossing_probability / poisson_pmf(n, n);
}
