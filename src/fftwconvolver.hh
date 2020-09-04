#ifndef __fftwconvolver_hh__
#define __fftwconvolver_hh__

#include <vector>
#include <complex>
#include <fftw3.h>

class FFTWConvolver {
public:
    FFTWConvolver(int maximum_input_size);
    ~FFTWConvolver();
    void convolve_same_size(int size, const double* input_a, const double* input_b, double* output);
private:
    int maximum_input_size;

    std::complex<double>* tmp_complex;


    // The r2c plans perform, for various sizes, a real to complex FFT with input at r2c_in and output at r2c_out
    double* r2c_in;
    std::complex<double>* r2c_out;
    std::vector<fftw_plan> r2c_plans;
    fftw_plan memoized_r2c_plan(int rounded_size);

    // The c2r plans perform, for various sizes, a complex to real FFT with input at c2r_in and output at c2r_out
    std::complex<double>* c2r_in;
    double* c2r_out;
    std::vector<fftw_plan> c2r_plans;
    fftw_plan memoized_c2r_plan(int rounded_size);

};

#endif
