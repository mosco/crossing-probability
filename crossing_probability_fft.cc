#include <cassert>
#include <cstring>
#include <vector>
extern "C" {
    #include <quadmath.h>
}
#include <complex.h>
#include <fftw3.h>

#include "fftw_wrappers.hh"

template<class T>
class FFTWConvolver {
public:
    FFTWConvolver(int maximum_input_size);
    ~FFTWConvolver();
    void convolve_same_size(int input_buffers_size, const T* input_buffer_0, const T* input_buffer_1, T* output_buffer);
private:
    int maximum_input_size;
    vector<FFTW_R2C_1D_Executor<T>*> r2c_executors;
    vector<FFTW_C2R_1D_Executor<T>*> c2r_executors;
    vector<typename fftw_traits<T>::complex_type> tmp;
};

template<class T>
FFTWConvolver<T>::FFTWConvolver(int maximum_input_size) :
    maximum_input_size(maximum_input_size),
    r2c_executors(maximum_input_size+1, NULL),
    c2r_executors(maximum_input_size+1, NULL),
    tmp(maximum_input_size+1)
{
}

template<class T>
void FFTWConvolver<T>::convolve_same_size(int input_buffers_size, const T* input_buffer_0, const T* input_buffer_1, T* output_buffer)
{
    int padded_length = 2*input_buffers_size;
    assert(input_buffers_size <= maximum_input_size);
    if (r2c_executors[input_buffers_size] == NULL) {
        //cerr << "Allocating " << input_buffers_size << endl;
        assert(r2c_executors[input_buffers_size] == NULL);
        assert(c2r_executors[input_buffers_size] == NULL);
        r2c_executors[input_buffers_size] = new FFTW_R2C_1D_Executor<T>(2*input_buffers_size);
        c2r_executors[input_buffers_size] = new FFTW_C2R_1D_Executor<T>(2*input_buffers_size);
    }
    //cerr << "accessing " << input_buffers_size << endl;

    FFTW_R2C_1D_Executor<T>& fft = *r2c_executors[input_buffers_size];
    assert (fft.input_size == padded_length);

    FFTW_C2R_1D_Executor<T>& ifft = *c2r_executors[input_buffers_size];
    assert (ifft.input_size == padded_length/2+1);

    fft.set_input_zeropadded(input_buffer_0, input_buffers_size);
    fft.execute();
    assert(static_cast<int>(tmp.size()) >= fft.output_size);
    memcpy(&tmp[0], fft.output_buffer, fft.output_size * sizeof(typename fftw_traits<T>::complex_type));

    fft.set_input_zeropadded(input_buffer_1, input_buffers_size);
    fft.execute();

    // Perform element-wise product of FFT(a) and FFT(b)
    // then compute inverse fourier transform.
    assert(fft.output_size == ifft.input_size);

    for (int i = 0; i < fft.output_size; ++i) {
        // FFTW returns unnormalized output. To normalize it one must divide each element
        // of the result by the number of elements.
        ifft.input_buffer[i] = tmp[i] * fft.output_buffer[i] / padded_length;
    }

    ifft.execute();

    assert(ifft.output_size == padded_length);
    memcpy(output_buffer, ifft.output_buffer, input_buffers_size * sizeof(T));
}

template<class T>
FFTWConvolver<T>::~FFTWConvolver()
{
    for (size_t i = 0; i < r2c_executors.size(); ++i) {
        if (r2c_executors[i] != NULL) {
            delete r2c_executors[i];
        }
        if (c2r_executors[i] != NULL) {
            delete c2r_executors[i];
        }
    }
}


template<class T>
T compute_crossing_probability_fft(const vector<double>& lower_bounds, const vector<double>& upper_bounds)
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

    vector<T> Qs0(n+1, 0.0);
    Qs0[0] = 1.0;
    vector<T> Qs1(n+1, 0.0);

    vector<T>* buffers[] = {&Qs0, &Qs1};

    T prev_location = 0.0;
    int lower_bound_count = 0;
    int upper_bound_count = 0;

    //vector<typename fftw_traits<T>::complex_type> tmp(2*n);
    FFTWConvolver<T> convolver(n+1);
    for (unsigned int i = 0; i < bounds.size(); ++i) {
        // cout << "--------------------------------------------\n";
        // cout << "Iteration " << i << endl;

        const vector<T>& from = *buffers[i % 2];
        vector<T>& to = *buffers[(i+1) % 2];

        //cout << "From:\n";
        //print_array(&from[0], n+1);
        //cout << endl;

        //cout << "lower_bound_count: " << lower_bound_count << endl;
        //cout << "upper_bound_count: " << upper_bound_count << endl;

        T location = bounds[i].location;

        int cur_size = lower_bound_count - upper_bound_count + 1;
        vector<T> pmf(cur_size);
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

    vector<T>& last_to_buffer = *buffers[bounds.size() % 2];
    T poisson_process_noncrossing_probability = last_to_buffer[n];
    T binomial_process_noncrossing_probability = poisson_process_noncrossing_probability / poisson_pmf(n, n);

    return 1-binomial_process_noncrossing_probability;
}
