#include <iostream>
#include <cassert>
#include <cstring>
#include "fftw_wrappers.hh"
#include "fftwconvolver.hh"
#include "aligned_mem.hh"

using namespace std;


// Constructing FFTW_R2C_1D_Executor objects for every possible input size is expensive,
// so we round up the sizes and construct an object for every N-th size.
const int ROUNDING = 128;
const int MINIMUM_SIZE_FOR_FFTW_CONVOLUTION = 80;

FFTWConvolver::FFTWConvolver(int maximum_input_size) :
    maximum_input_size(maximum_input_size+ROUNDING-1),
    r2c_executors(maximum_input_size+ROUNDING, NULL),
    c2r_executors(maximum_input_size+ROUNDING, NULL)
{
    tmp_complex = allocate_aligned_complexes(maximum_input_size);
    tmp_double_0 = allocate_aligned_doubles(maximum_input_size);
    tmp_double_1 = allocate_aligned_doubles(maximum_input_size);
    tmp_double_2 = allocate_aligned_doubles(maximum_input_size);
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

void elementwise_complex_product(int size, const complex<double>* __restrict__ src0, const complex<double>* __restrict__  src1, complex<double>* __restrict__ dest, double multiplicative_constant)
{
    for (int i = 0; i < size; ++i) {
        dest[i] = multiplicative_constant*src0[i]*src1[i];
    }
}

void FFTWConvolver::convolve_same_size(int input_buffers_size, const double* input_buffer_0, const double* input_buffer_1, double* output_buffer)
{
    assert(input_buffers_size <= maximum_input_size);
    if (input_buffers_size < MINIMUM_SIZE_FOR_FFTW_CONVOLUTION) {
        memcpy(tmp_double_0, input_buffer_0, input_buffers_size*sizeof(double));
        memcpy(tmp_double_1, input_buffer_1, input_buffers_size*sizeof(double));
        convolve_same_size_naive(input_buffers_size, tmp_double_0, tmp_double_1, tmp_double_2);
        memcpy(output_buffer, tmp_double_2, input_buffers_size*sizeof(double));
        return;
    }
    int padded_length = ((2*input_buffers_size+ROUNDING-1)/ROUNDING)*ROUNDING;
    int executor_index = padded_length/ROUNDING - 1;

    assert(input_buffers_size <= maximum_input_size);
    if (r2c_executors[executor_index] == NULL) {
        assert(r2c_executors[executor_index] == NULL);
        assert(c2r_executors[executor_index] == NULL);
        r2c_executors[executor_index] = new FFTW_R2C_1D_Executor(padded_length);
        c2r_executors[executor_index] = new FFTW_C2R_1D_Executor(padded_length);
    }

    FFTW_R2C_1D_Executor& fft = *r2c_executors[executor_index];
    assert (fft.input_size >= padded_length);

    FFTW_C2R_1D_Executor& ifft = *c2r_executors[executor_index];
    assert (ifft.input_size >= padded_length/2+1);

    fft.set_input_zeropadded(input_buffer_0, input_buffers_size);
    fft.execute();
    std::memcpy(tmp_complex, fft.output_buffer, fft.output_size * sizeof(complex<double>));

    fft.set_input_zeropadded(input_buffer_1, input_buffers_size);
    fft.execute();

    // Perform element-wise product of FFT(a) and FFT(b)
    // then compute inverse fourier transform.
    assert(fft.output_size == ifft.input_size);
    // FFTW returns unnormalized output. To normalize it one must divide each element
    // of the result by the number of elements.
    elementwise_complex_product(fft.output_size, tmp_complex, fft.output_buffer, ifft.input_buffer, 1.0/double(padded_length));

    ifft.execute();

    assert(ifft.output_size == padded_length);
    std::memcpy(output_buffer, ifft.output_buffer, input_buffers_size * sizeof(double));
}

FFTWConvolver::~FFTWConvolver()
{
    for (size_t i = 0; i < r2c_executors.size(); ++i) {
        if (r2c_executors[i] != NULL) {
            delete r2c_executors[i];
        }
        if (c2r_executors[i] != NULL) {
            delete c2r_executors[i];
        }
    }
    free_aligned_mem(tmp_double_0);
    free_aligned_mem(tmp_double_1);
    free_aligned_mem(tmp_double_2);
    free_aligned_mem(tmp_complex);
}
