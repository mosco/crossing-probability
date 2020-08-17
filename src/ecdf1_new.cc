#include <vector>
#include <iostream>
#include <stdexcept>

#include "ecdf1_new.hh"
#include "common.hh"
#include "poisson_pmf.hh"
#include "fftwconvolver.hh"
#include "aligned_mem.hh"
#include "string_utils.hh"

using namespace std;

vector<double> poisson_B_noncrossing_probability_n2(int n, double intensity, const vector<double>& B, int jump_size)
{
    assert(jump_size <= n);
    DoubleBuffer<double> buffers(n+1, 0.0);
    DoubleBuffer<double> minibuffers(jump_size, 0.0);
    buffers.get_src()[0] = 1.0;

    FFTWConvolver fftconvolver(n+1);

    double* pmf = allocate_aligned_doubles(n+1);
    PoissonPMFGenerator pmfgen(n+1);

    double* tmp = allocate_aligned_doubles(n+1);

    int n_steps = B.size();
    double I_prev_location = 0.0;
    int I_prev = -1;
    int I = I_prev+jump_size;
    while (true) {
        // cout << "I: " << I << " I_prev: " << I_prev << endl;
        pmfgen.compute_array(n-I_prev, intensity*(B[I]-I_prev_location));
        fftconvolver.convolve_same_size(n-I_prev, pmfgen.get_array(), &buffers.get_src()[I_prev+1], tmp);
        fill(&buffers.get_dest()[0], &buffers.get_dest()[I+1], 0.0);
        copy(&tmp[I-I_prev], &tmp[n-I_prev], &buffers.get_dest()[I+1]);

        copy(&buffers.get_src()[I_prev+1], &buffers.get_src()[I+1], &minibuffers.get_src()[0]);
        double i_prev_location = I_prev_location;
        for (int i = I_prev+1; i < I; i++) {
            pmfgen.compute_array(I-i+1, intensity*(B[i]-i_prev_location));
            fftconvolver.convolve_same_size(I-i+1, pmfgen.get_array(), &minibuffers.get_src()[i-I_prev-1], &minibuffers.get_dest()[i-I_prev-1]);

            double prob_exit_now = minibuffers.get_dest()[i-I_prev-1];
            double lambda = intensity*(B[I]- B[i]);
            for (int j = I+1; j < n+1; ++j) {
                buffers.get_dest()[j] -= prob_exit_now * pmfgen.evaluate_pmf(lambda, j-i);
            }

            minibuffers.get_dest()[i-I_prev-1] = 0.0;
            minibuffers.get_src()[i-I_prev-1] = 0.0;
            minibuffers.flip();
            i_prev_location = B[i];
        }
        // cout << "I: " << I << " I_prev: " << I_prev << endl;
        I_prev = I;
        // cout << "I: " << I << " I_prev: " << I_prev << endl;
        I_prev_location = B[I];
        buffers.flip();
        if (I == (n_steps-1)) {
            break;
        }

        // cout << "I: " << I << " I_prev: " << I_prev << endl;
        I = min(I+jump_size, n_steps-1);
        // cout << "I: " << I << " I_prev: " << I_prev << endl;
    }
    pmfgen.compute_array(n-n_steps+1, intensity*(1.0-I_prev_location));
    fftconvolver.convolve_same_size(n-n_steps+1, pmfgen.get_array(), &buffers.get_src()[n_steps], &buffers.get_dest()[n_steps]);
    fill(&buffers.get_dest()[0], &buffers.get_dest()[n_steps], 0.0);

    free_aligned_mem(tmp);
    free_aligned_mem(pmf);
    return buffers.get_dest();
}

double ecdf1_new_B(int n, const vector<double>& B)
{
    if ((int)B.size() != n) {
        stringstream ss;
        ss << "Expecting exactly " << n << "bounds. " << "Got " << B.size() << ".";
        throw runtime_error(ss.str() + vector_to_string(B));
    }
    // Asymptotically any k in the range [logn, n/logn] should give optimal results as n goes to infinity.
    // Setting k=c*sqrt(n) and minimizing the asymptotic runtime, we obtain k=sqrt(2*n),
    // however, empirically 5*sqrt(n) gives better results.
    int k = min(5*sqrt(n),n/log2(n));
    vector<double> poisson_nocross_probabilities = poisson_B_noncrossing_probability_n2(n, n, B, k);
    return poisson_nocross_probabilities[n] / poisson_pmf(n, n);
}
// For n=10000, best results k=400...600

double ecdf1_new_b(int n, const vector<double>& b)
{
    if ((int)b.size() != n) {
        stringstream ss;
        ss << "Expecting exactly " << n << "bounds. " << "Got " << b.size() << ".";
        throw runtime_error(ss.str());
    }

    vector<double> symmetric_steps(n, 0.0);
    for (int i = n-b.size(); i < n; ++i) {
        symmetric_steps[i] = 1.0 - b[b.size() - 1 - i];
    }

    return ecdf1_new_B(n, symmetric_steps);
}
