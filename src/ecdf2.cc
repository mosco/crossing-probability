#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <ctime>

#include "ecdf2.hh"
#include "fftwconvolver.hh"
#include "aligned_mem.hh"
#include "common.hh"
#include "poisson_pmf.hh"
#include "string_utils.hh"

using namespace std;

enum BoundType {bSTEP, BSTEP, END};  

struct Bound {
    double location;
    BoundType tag;
};

// Needed for using sort()
static bool operator<(Bound b0, Bound b1)
{
    return (b0.location < b1.location);
}

static vector<Bound> join_all_bounds(const vector<double>& b, const vector<double>& B)
{
    vector<Bound> bounds;
    bounds.reserve(b.size()+B.size()+1);

    Bound bound;

    for (int i = 0; i < (int)b.size(); ++i) {
        bound.location = b[i];
        bound.tag = bSTEP;
        bounds.push_back(bound);
    }

    for (int i = 0; i < (int)B.size(); ++i) {
        bound.location = B[i];
        bound.tag = BSTEP;
        bounds.push_back(bound);
    }

    sort(bounds.begin(), bounds.end());

    bound.location = 1.0;
    bound.tag = END;
    bounds.push_back(bound);

    return bounds;
}

void update_buffers_and_step_counts(Bound bound, DoubleBuffer<double>& buffers, int& b_step_count, int& B_step_count)
{
    if (bound.tag == bSTEP) {
        ++b_step_count;
        buffers.get_dest()[b_step_count] = 0.0;
        // The following line is not necessary, but keeps the arrays cleaner.
        // This is helpful for doing debug prints.
        //buffers.get_src()[b_step_count] = 0.0; 
    } else if (bound.tag == BSTEP) {
        buffers.get_dest()[B_step_count] = 0.0;
        // The following line is not necessary, but keeps the arrays cleaner.
        // This is helpful for doing debug prints.
        ////buffers.get_src()[B_step_count] = 0.0;
        ++B_step_count;
    } else {
        if (bound.tag != END) {
            cout << "tag: " << bound.tag << "\n";
            throw runtime_error("Expecting END tag");
        }
    }
}

// TODO: Split function into 2 cases: with_fft and no_fft
vector<double> poisson_process_noncrossing_probability(int n, double intensity, const vector<double>& b, const vector<double>& B, bool use_fft)
{
    vector<Bound> bounds = join_all_bounds(b, B);

    DoubleBuffer<double> buffers(n+1, 0.0);
    buffers.get_src()[0] = 1.0;

    FFTWConvolver fftconvolver(n+1); 

    PoissonPMFGenerator pmfgen(n+1);

    int b_step_count = 0;
    int B_step_count = 0;

    double prev_location = 0.0;

    for (unsigned int i = 0; i < bounds.size(); ++i) {
        int cur_size = b_step_count - B_step_count + 1;

        double lambda = intensity*(bounds[i].location-prev_location);
        pmfgen.compute_array(cur_size, lambda);

        if (use_fft) {
            fftconvolver.convolve_same_size(cur_size, pmfgen.get_array(), &buffers.get_src()[B_step_count], &buffers.get_dest()[B_step_count]);
        } else {
            convolve_same_size(cur_size, pmfgen.get_array(), &buffers.get_src()[B_step_count], &buffers.get_dest()[B_step_count]);
        }

        update_buffers_and_step_counts(bounds[i], buffers, b_step_count, B_step_count);

        prev_location = bounds[i].location;
        buffers.flip();
    }

    return buffers.get_src();
}

double ecdf2(int n, const vector<double>& b, const vector<double>& B, bool use_fft)
{
    //TODO: Check that the boundaries are valid here!
    //Python extensions call this function directly.

    if (int(b.size()) != n) {
        throw runtime_error("Expecting exactly n bounds: b_1, ..., b_n");
    }
    if (int(B.size()) != n) {
        throw runtime_error("Expecting exactly n bounds: B_1, ..., B_n");
    }

    vector<double> poisson_nocross_probs = poisson_process_noncrossing_probability(n, n, b, B, use_fft);

    return poisson_nocross_probs[n] / poisson_pmf(n, n);
}

