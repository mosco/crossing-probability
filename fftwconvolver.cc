#include <cassert>
#include <cstring>
#include "fftw_wrappers.hh"
#include "fftwconvolver.hh"

// Constructing FFTW_R2C_1D_Executor objects for every input size takes a lot of time, so we round up the sizes.
const int ROUNDING = 128;

FFTWConvolver::FFTWConvolver(int maximum_input_size) :
    maximum_input_size(maximum_input_size+ROUNDING-1),
    r2c_executors(maximum_input_size+ROUNDING, NULL),
    c2r_executors(maximum_input_size+ROUNDING, NULL),
    tmp(maximum_input_size+ROUNDING)
{
}

void FFTWConvolver::convolve_same_size(int input_buffers_size, const double* input_buffer_0, const double* input_buffer_1, double* output_buffer)
{
    int size_rounded = ROUNDING*((input_buffers_size+ROUNDING-1)/ROUNDING);
    int padded_length = 2*size_rounded;
    assert(size_rounded <= maximum_input_size);
    if (r2c_executors[size_rounded/ROUNDING] == NULL) {
        //cerr << "Allocating " << size_rounded << endl;
        assert(r2c_executors[size_rounded/ROUNDING] == NULL);
        assert(c2r_executors[size_rounded/ROUNDING] == NULL);
        r2c_executors[size_rounded/ROUNDING] = new FFTW_R2C_1D_Executor(padded_length);
        c2r_executors[size_rounded/ROUNDING] = new FFTW_C2R_1D_Executor(padded_length);
    }
    //cerr << "accessing " << size_rounded << endl;

    FFTW_R2C_1D_Executor& fft = *r2c_executors[size_rounded/ROUNDING];
    assert (fft.input_size == padded_length);

    FFTW_C2R_1D_Executor& ifft = *c2r_executors[size_rounded/ROUNDING];
    assert (ifft.input_size == padded_length/2+1);

    fft.set_input_zeropadded(input_buffer_0, input_buffers_size);
    fft.execute();
    assert(static_cast<int>(tmp.size()) >= fft.output_size);
    std::memcpy(&tmp[0], fft.output_buffer, fft.output_size * sizeof(double complex));

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
}
