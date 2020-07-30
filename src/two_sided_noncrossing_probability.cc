#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <ctime>
#include "fftwconvolver.hh"
#include "aligned_mem.hh"
#include "common.hh"
#include "poisson_pmf.hh"
#include "string_utils.hh"

using namespace std;

enum BoundType {START, UBSTEP, LBSTEP, END};

struct Bound {
    double location;
    BoundType tag;
    int lb, ub; // Height of lower boundary and upper boundary functions immediately before the current step
};

static bool operator<(Bound b0, Bound b1)
{
    return (b0.location < b1.location);
}

static vector<Bound> join_all_bounds(const vector<double>& ubsteps, const vector<double>& lbsteps)
{
    assert(ubsteps.size() >= lbsteps.size());

    vector<Bound> bounds;
    bounds.reserve(ubsteps.size()+lbsteps.size()+2);
    Bound b;

    b.location = 0.0;
    b.tag = START;
    bounds.push_back(b);

    for (int i = 0; i < (int)ubsteps.size(); ++i) {
        b.location = ubsteps[i];
        b.tag = UBSTEP;
        bounds.push_back(b);
    }

    for (int i = 0; i < (int)lbsteps.size(); ++i) {
        b.location = lbsteps[i];
        b.tag = LBSTEP;
        bounds.push_back(b);
    }

    sort(bounds.begin()+1, bounds.end());

    b.location = 1.0;
    b.tag = END;
    bounds.push_back(b);

    int lb = -1;
    int ub = 1;
    for (int i = 0; i < bounds.size(); ++i) {
        bounds[i].lb = lb;
        bounds[i].ub = ub;
        if (bounds[i].tag == UBSTEP) {
            ++ub;
        } else if (bounds[i].tag == LBSTEP) {
            ++lb;
        }
    }

    return bounds;
}

static bool lower_and_upper_boundaries_cross(const vector<double>& lbsteps, const vector<double>& ubsteps)
{
    if (lbsteps.size() > ubsteps.size()) {
        cout << "The lower and upper boundaries cross: lower(1) > upper(1).\n";
        return true;
    }
    for (size_t i = 0; i < lbsteps.size(); ++i) {
        if (lbsteps[i] < ubsteps[i]) {
            cout << "The lower and upper boundaries cross! i=" << i << ".\n";
            return true;
        }
    }
    return false;
}

vector<double> poisson_process_noncrossing_probability(double intensity, const vector<double>& g_steps, const vector<double>& h_steps, bool use_fft)
{
    if (lower_and_upper_boundaries_cross(g_steps, h_steps)) {
        return vector<double>();
    }

    vector<Bound> bounds = join_all_bounds(h_steps, g_steps);
    int n = h_steps.size();
    DoubleBuffer<double> buffers(n+1, 0.0);
    buffers.get_src()[0] = 1.0;
    FFTWConvolver fftconvolver(n+1);
    PoissonPMFGenerator pmfgen(n+1);
    int h_step_count = 0;
    int g_step_count = 0;
    double prev_location = 0.0;

    for (unsigned int i = 0; i < bounds.size(); ++i) {
        int cur_size = h_step_count - g_step_count + 1;

        double lambda = intensity*(bounds[i].location-prev_location);
        pmfgen.compute_array(cur_size, lambda);

        if (use_fft) {
            fftconvolver.convolve_same_size(cur_size, pmfgen.get_array(), &buffers.get_src()[g_step_count], &buffers.get_dest()[g_step_count]);
        } else {
            convolve_same_size(cur_size, pmfgen.get_array(), &buffers.get_src()[g_step_count], &buffers.get_dest()[g_step_count]);
        }

        BoundType tag = bounds[i].tag;
        if (tag == UBSTEP) {
            ++h_step_count;
            buffers.get_dest()[h_step_count] = 0.0;
            buffers.get_src()[h_step_count] = 0.0; // Not strictly necessary. This just keeps the arrays cleaner when printing.
        } else if (tag == LBSTEP) {
            buffers.get_dest()[g_step_count] = 0.0;
            buffers.get_src()[g_step_count] = 0.0; // Not strictly necessary. This just keeps the arrays cleaner when printing.
            ++g_step_count;
        } else {
            assert(tag == END);
            break;
        }
        prev_location = bounds[i].location;
        buffers.flip();
    }

    return buffers.get_dest();
}

double ecdf_noncrossing_probability(int n, const vector<double>& lowerbound_steps, const vector<double>& upperbound_steps, bool use_fft)
{
    if ((int)lowerbound_steps.size() > n) {
        stringstream ss;
        ss << "Empirical CDF must cross lower boundary g(t) since g(1)==" << lowerbound_steps.size() << " > n and the number of samples is n." << endl;
        throw runtime_error(ss.str());
    }
    vector<double> processed_upperbound_steps(n, 0.0);
    if (upperbound_steps.size() == 0) {
        // Special case, only the lower bound is specified.
        // We treat this as an implicit upper bound satisfying h(t) = n for all t.
    } else {
        if (lower_and_upper_boundaries_cross(lowerbound_steps, upperbound_steps)) {
            return 0.0;
        }
        if ((int)upperbound_steps.size() < n) {
            stringstream ss;
            ss << "Empirical CDF must cross lower boundary g(t) since h(1)==" << upperbound_steps.size() << " > n and the number of samples is n. h_steps:" << endl;
            throw runtime_error(ss.str() + vector_to_string(upperbound_steps));
        }
        copy(upperbound_steps.begin(), upperbound_steps.begin() + n, processed_upperbound_steps.begin());
    }

    clock_t start_clock = clock();
    vector<double> poisson_nocross_probs = poisson_process_noncrossing_probability(n, lowerbound_steps, processed_upperbound_steps, use_fft);
    clock_t end_clock = clock();

    cout << "Elapsed clock: " << (end_clock-start_clock)/double(CLOCKS_PER_SEC) << endl;
    return poisson_nocross_probs[n] / poisson_pmf(n, n);
}

//void propagate_poisson_process(DoubleBuffer<double>& buffers, const vector<Bound> bounds, int start, int end, FFTWConvolver& fftconvolver, PoissonPMFGenerator& pmfgen)
//{
//    assert(0 <= start);
//    assert(start < end);
//    assert(end <= bounds.size());
//
//    for (unsigned int i = start; i < end; ++i) {
//        // Propagate poisson probabilities from Q_i to Q_{i+1}
//        // The input range is (lb+1, ..., ub-1) inclusive
//        // If the current step is an upper boundary step, the output range is (lb+1, ..., ub)
//        // If it is a lower boundary step, the output range is (lb+2, ..., ub-1)
//        // In both cases we compute the output for (lb+1, ..., ub) and then truncate.
//        int cur_size = bounds[i].ub - bounds[i].lb;
//
//        double lambda = intensity*(bounds[i+1].location - bounds[i].location);
//        pmfgen.compute_array(cur_size, lambda);
//
//        fftconvolver.convolve_same_size(cur_size, pmfgen.get_array(), &buffers.get_src()[bounds[i].lb+1], &buffers.get_dest()[bounds[i].lb+1]);
//
//        //************ TODO: CONTINUE ****************
//
//        if (bounds[i].tag == UBSTEP) {
//            ++h_step_count;
//            buffers.get_dest()[h_step_count] = 0.0;
//            buffers.get_src()[h_step_count] = 0.0; // Not strictly necessary. This just keeps the arrays cleaner when printing.
//        } else if (bounds[i].tag == G_STEP) {
//            buffers.get_dest()[g_step_count] = 0.0;
//            buffers.get_src()[g_step_count] = 0.0; // Not strictly necessary. This just keeps the arrays cleaner when printing.
//            ++g_step_count;
//        } else {
//            assert(tag == END);
//            break;
//        }
//        buffers.flip();
//    }
//
//    return buffers.get_dest();
//}
//
//void advance_multistep(DoubleBuffer<double>& buffers, const vector<Bound> bounds, int start, int end, FFTWConvolver& fftconvolver, PoissonPMFGenerator& pmfgen, double* pmf)
//{
//    if (bounds[end].g_before >= bounds[start].h_before) {
//        propagate_poisson_process(buffers, bounds, start, end, fftconvolver, pmfgen, );
//    } else {
//        TODO: Fancy computation
//    }
//}
//
//vector<double> poisson_process_noncrossing_probability_n2(double intensity, const vector<double>& g_steps, const vector<double>& h_steps, int multistep_size)
//{
//    if (lower_and_upper_boundaries_cross(g_steps, h_steps)) {
//        return vector<double>();
//    }
//
//    vector<Bound> bounds = join_all_bounds(h_steps, g_steps);
//    int n = h_steps.size();
//    DoubleBuffer<double> buffers(n+1, 0.0);
//    buffers.get_src()[0] = 1.0;
//
//    FFTWConvolver fftconvolver(n+1);
//    PoissonPMFGenerator pmfgen(n+1);
//    double* pmf = allocate_aligned_doubles(n+1);
//
//    int multistep_start = 1;
//    while (true) {
//        int multistep_end = min(multistep_start + multistep_size, bounds.size());
//        advance_multistep(buffers, bounds, multistep_start, multistep_end, fftconvolver, pmfgen, pmf);
//        if (multistep_end == bounds.size()) {
//            break;
//        }
//        buffers.flip()
//        multistep_start = multistep_end;
//    }
//    free_aligned_mem(pmf);
//    return buffers.get_dest();
//}
//
//double ecdf_noncrossing_probability_n2(int n, const vector<double>& g_steps, const vector<double>& h_steps, bool use_fft)
//{
//    if ((int)g_steps.size() > n) {
//        stringstream ss;
//        ss << "Empirical CDF must cross lower boundary g(t) since g(1)==" << g_steps.size() << " > n and the number of samples is n." << endl;
//        throw runtime_error(ss.str());
//    }
//    vector<double> processed_h_steps(n, 0.0);
//    if (h_steps.size() == 0) {
//        // Special case, only the lower bound is specified.
//        // We treat this as an implicit upper bound satisfying h(t) = n for all t.
//    } else {
//        if (lower_and_upper_boundaries_cross(g_steps, h_steps)) {
//            return 0.0;
//        }
//        if ((int)h_steps.size() < n) {
//            stringstream ss;
//            ss << "Empirical CDF must cross lower boundary g(t) since h(1)==" << h_steps.size() << " > n and the number of samples is n. h_steps:" << endl;
//            throw runtime_error(ss.str() + vector_to_string(h_steps));
//        }
//        copy(h_steps.begin(), h_steps.begin() + n, processed_h_steps.begin());
//    }
//
//    clock_t start_clock = clock();
//    vector<double> poisson_nocross_probs = poisson_process_noncrossing_probability_n2(n, g_steps, processed_h_steps, use_fft);
//    clock_t end_clock = clock();
//    cout << "Elapsed clock: " << (end_clock-start_clock)/double(CLOCKS_PER_SEC) << endl;
//    return poisson_nocross_probs[n] / poisson_pmf(n, n);
//}


