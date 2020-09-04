#include <iostream>
#include <cassert>
#include <cstring>
#include "fftwconvolver.hh"
#include "aligned_mem.hh"

using namespace std;


// const int ROUNDING = 128;
// const int MINIMUM_SIZE_FOR_FFTW_CONVOLUTION = 80;

// Good for ecdf2-mn2017 of one-sided boundary
//const int ROUNDING = 2048;
//const int MINIMUM_SIZE_FOR_FFTW_CONVOLUTION = 80;

// Good for ecdf1-new
const int ROUNDING = 2048;
const int MINIMUM_SIZE_FOR_FFTW_CONVOLUTION = 128;


int round_up(int n, int rounding)
{
    assert(rounding >= 0);
    return ((n+rounding-1)/rounding)*rounding;
}

FFTWConvolver::FFTWConvolver(int maximum_input_size) :
    maximum_input_size(maximum_input_size+ROUNDING-1),
    r2c_plans(round_up(2*maximum_input_size, ROUNDING)/ROUNDING, NULL),
    c2r_plans(round_up(2*maximum_input_size, ROUNDING)/ROUNDING, NULL)
{
    int maximum_padded_input_size = round_up(2*maximum_input_size, ROUNDING);

    r2c_in = allocate_aligned_doubles(maximum_padded_input_size);
    r2c_out = allocate_aligned_complexes(maximum_padded_input_size);

    c2r_in = allocate_aligned_complexes(maximum_padded_input_size);
    c2r_out = allocate_aligned_doubles(maximum_padded_input_size);

    tmp_complex = allocate_aligned_complexes(maximum_padded_input_size);
}

void convolve_same_size_naive(int size, const double* __restrict__ src0, const double* __restrict__ src1, double* __restrict__ dest)
{
    for (int j = 0; j < size; ++j) {
        double convolution_at_j = 0.0;
        for (int k = 0; k <= j; ++k) {
            convolution_at_j += src0[k] * src1[j-k];
        }
        dest[j] = convolution_at_j;
        
    }
}

void elementwise_complex_product(
    int size,
    const complex<double>* __restrict__ src0,
    const complex<double>* __restrict__  src1,
    complex<double>* __restrict__ dest,
    double multiplicative_constant)
{
    for (int i = 0; i < size; ++i) {
        dest[i] = multiplicative_constant*src0[i]*src1[i];
    }
}

fftw_plan FFTWConvolver::memoized_r2c_plan(int rounded_size)
{
    assert(rounded_size > 0);
    assert((rounded_size % ROUNDING) == 0);

    int index = rounded_size/ROUNDING - 1;
    assert(index < r2c_plans.size());

    if (r2c_plans[index] == NULL) {
        r2c_plans[index] = fftw_plan_dft_r2c_1d(rounded_size, r2c_in, reinterpret_cast<fftw_complex*>(r2c_out), FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
    }

    return r2c_plans[index];
}

fftw_plan FFTWConvolver::memoized_c2r_plan(int rounded_size)
{
    assert(rounded_size > 0);
    assert((rounded_size % ROUNDING) == 0);

    int index = rounded_size/ROUNDING - 1;
    assert(index < c2r_plans.size());

    if (c2r_plans[index] == NULL) {
        c2r_plans[index] = fftw_plan_dft_c2r_1d(rounded_size, reinterpret_cast<fftw_complex*>(c2r_in), c2r_out, FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
    }

    return c2r_plans[index];
}

template<class T>
void copy_zero_padded(const T* src, T* dest, int src_size, int dest_size)
{
    memcpy(dest, src, sizeof(T)*src_size);
    memset(&dest[src_size], 0, sizeof(T)*(dest_size-src_size));
}

void FFTWConvolver::convolve_same_size(
    int size,
    const double* __restrict__ input_a,
    const double* __restrict__ input_b,
    double* __restrict__ output)
{
    if (size > maximum_input_size) {
        stringstream ss;
        ss << "FFTWConvolver::convolve_same_size received input of size " << size << ". This is bigger than maximum_input_size==" << maximum_input_size;
        throw runtime_error(ss.str());
    }
    if (size <= 0) {
        return; // Nothing to do
    }

    if (size < MINIMUM_SIZE_FOR_FFTW_CONVOLUTION) {
        convolve_same_size_naive(size, input_a, input_b, output);
        return;
    }

    int padded_size = round_up(2*size, ROUNDING);
    
    // tmp_complex <- FFT(zeropad(input_a));
    copy_zero_padded(input_a, r2c_in, size, padded_size);
    //fftw_execute(memoized_r2c_plan(padded_size));
    //memcpy(tmp_complex, r2c_out, sizeof(complex<double>)*(padded_size/2 + 1));
    fftw_execute_dft_r2c(memoized_r2c_plan(padded_size), r2c_in, reinterpret_cast<fftw_complex*>(tmp_complex));
    // Try this instead of last two lines:
    // fftw_execute_dft_r2c(memoized_r2c_plan(padded_length), r2c_in, tmp_complex);

    // r2c_out <- FFT(zeropad(input_b));
    copy_zero_padded(input_b, r2c_in, size, padded_size); 
    fftw_execute(memoized_r2c_plan(padded_size));

    // Perform element-wise product of FFT(a) and FFT(b) and then compute inverse fourier transform.
    // FFTW returns unnormalized output. To normalize it one must divide each element of the result by the number of elements.
    elementwise_complex_product(padded_size/2 + 1, tmp_complex, r2c_out, c2r_in, 1.0/double(padded_size));
    fftw_execute(memoized_c2r_plan(padded_size));
    std::memcpy(output, c2r_out, size * sizeof(double));
}

FFTWConvolver::~FFTWConvolver()
{
    for (size_t i = 0; i < r2c_plans.size(); ++i) {
        if (r2c_plans[i] != NULL) {
            fftw_destroy_plan(r2c_plans[i]);
        }
    }

    for (size_t i = 0; i < c2r_plans.size(); ++i) {
        if (c2r_plans[i] != NULL) {
            fftw_destroy_plan(c2r_plans[i]);
        }
    }

    free_aligned_mem(r2c_in);
    free_aligned_mem(r2c_out);

    free_aligned_mem(c2r_in);
    free_aligned_mem(c2r_out);

    free_aligned_mem(tmp_complex);
}
