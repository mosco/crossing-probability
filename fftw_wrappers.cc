
#include <cassert>
#include <cstring>
#include "fftw_wrappers.hh"

using namespace std;

FFTW_R2C_1D_Executor::FFTW_R2C_1D_Executor(int n_real_samples) :
    input_size(n_real_samples),
    input_buffer(fftwq_alloc_real(n_real_samples)),
    output_size(n_real_samples/2 + 1),
    output_buffer(fftwq_alloc_complex(n_real_samples/2 + 1))
{
    plan = fftwq_plan_dft_r2c_1d(n_real_samples, input_buffer, output_buffer, FFTW_ESTIMATE);
}

FFTW_R2C_1D_Executor::~FFTW_R2C_1D_Executor()
{
    fftwq_destroy_plan(plan);
    fftwq_free(input_buffer);
    fftwq_free(output_buffer);
}

void FFTW_R2C_1D_Executor::set_input_zeropadded(const __float128* buffer, int size)
{
    assert(size <= input_size);
    memcpy(input_buffer, buffer, sizeof(__float128)*size);
    memset(&input_buffer[size], 0, sizeof(__float128)*(input_size - size));
}

void FFTW_R2C_1D_Executor::set_input_zeropadded(const vector<__float128>& vec)
{
    set_input_zeropadded(&vec[0], vec.size());
}

void FFTW_R2C_1D_Executor::execute()
{
    fftwq_execute(plan);
}

vector<__float128 complex> FFTW_R2C_1D_Executor::get_output()
{
    return vector<__float128 complex>(output_buffer, output_buffer + output_size);
}

FFTW_C2R_1D_Executor::FFTW_C2R_1D_Executor(int n_real_samples) : 
    input_size(n_real_samples/2 + 1),
    input_buffer(fftwq_alloc_complex(n_real_samples/2 + 1)),
    output_size(n_real_samples),
    output_buffer(fftwq_alloc_real(n_real_samples))
{
    plan = fftwq_plan_dft_c2r_1d(n_real_samples, input_buffer, output_buffer, FFTW_ESTIMATE);
}

FFTW_C2R_1D_Executor::~FFTW_C2R_1D_Executor()
{
    fftwq_destroy_plan(plan);
    fftwq_free(input_buffer);
    fftwq_free(output_buffer);
}

void FFTW_C2R_1D_Executor::set_input(const __float128 complex* buffer, int size)
{
    assert(size == input_size);
    memcpy(input_buffer, buffer, sizeof(__float128 complex)*size);
    memset(&input_buffer[size], 0, sizeof(__float128 complex)*(input_size - size));
}

void FFTW_C2R_1D_Executor::set_input(const vector<__float128 complex>& vec)
{
    set_input(&vec[0], vec.size());
}

void FFTW_C2R_1D_Executor::execute()
{
    fftwq_execute(plan);
}

vector<__float128> FFTW_C2R_1D_Executor::get_output()
{
    return vector<__float128>(output_buffer, output_buffer + output_size);
}
