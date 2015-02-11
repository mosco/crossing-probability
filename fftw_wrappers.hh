#ifndef __FFTWWrapper_h__
#define __FFTWWrapper_h__

// Typical use of FFTW entails the following steps:
// 1. Allocate input and output buffers.
// 2. Compute a "plan" struct that. This tells FFTW what algorithms it should use when actually computing the FFT.
// 3. Execute the FFT/IFFT operation on the buffers from step 1.
//
// This file contains two classes that wrap this process nicely.
// Currently only one-dimensional transforms of real data to complex and back are supported.

#include <vector>
#include <complex.h>
#include <fftw3.h>
#include <cassert>

#include <iostream>

// Usage: (after initializing the class)
// 1. Fill input_buffer with input containing n_real_samples __float128 numbers
//    (note, set_input_zeropadded will copy your buffer with optional zero padding)
// 2. Run execute().
// 3. Extract output by calling get_output() or directly access output_buffer[0], ..., output_buffer[output_size-1].
//    Note that the output is composed of n_real_samples/2 + 1 complex numbers.
// 
// These 3 steps can be repeated many times.

template <class T>
struct fftw_traits {};

template<>
struct fftw_traits<double> {
    typedef fftw_plan plan_type;
    typedef double complex complex_type;
    
    static double* alloc_real(size_t n) 
    {
        return fftw_alloc_real(n); 
    }
    static complex_type* alloc_complex(size_t n) 
    {
        return fftw_alloc_complex(n); 
    }
    static fftw_plan plan_dft_r2c_1d(int n0, double* in, double complex* out, unsigned flags) 
    { 
        return fftw_plan_dft_r2c_1d(n0, in, out, flags); 
    }
    static fftw_plan plan_dft_c2r_1d(int n0, double complex* in, double* out, unsigned flags) 
    { 
        return fftw_plan_dft_c2r_1d(n0, in, out, flags); 
    }
    static void destroy_plan(fftw_plan plan) 
    { 
        fftw_destroy_plan(plan); 
    }
    static void free(void* p)
    {
        fftw_free(p);
    }
    static void execute(const fftw_plan plan)
    {
        fftw_execute(plan);
    }
};

template<>
struct fftw_traits<__float128> {
    typedef fftwq_plan plan_type;
    typedef __float128 complex complex_type;

    static __float128* alloc_real(size_t n)
    {
        return fftwq_alloc_real(n);
    }
    static __float128 complex* alloc_complex(size_t n) 
    {
        return fftwq_alloc_complex(n); 
    }
    static plan_type plan_dft_r2c_1d(int n0, __float128* in, __float128 complex* out, unsigned flags)
    {
        return fftwq_plan_dft_r2c_1d(n0, in, out, flags);
    }
    static fftwq_plan plan_dft_c2r_1d(int n0, __float128 complex* in, __float128* out, unsigned flags) 
    { 
        return fftwq_plan_dft_c2r_1d(n0, in, out, flags); 
    }
    static void destroy_plan(fftwq_plan plan)
    {
        fftwq_destroy_plan(plan); 
    }
    static void free(void* p)
    {
        fftwq_free(p);
    }
    static void execute(const fftwq_plan plan)
    {
        fftwq_execute(plan);
    }
};

template <class T>
class FFTW_R2C_1D_Executor {
public:
    typedef fftw_traits<T> traits;

    FFTW_R2C_1D_Executor(int n_real_samples);
    ~FFTW_R2C_1D_Executor();
    void set_input_zeropadded(const T* buffer, int size);
    void set_input_zeropadded(const std::vector<T>& vec);
    void execute();
    std::vector<typename traits::complex_type> get_output();

    const int input_size;
    T* const input_buffer;

    const int output_size;
    typename traits::complex_type* const output_buffer;

private:
    typename traits::plan_type plan;
};

// Usage of this class is similar to that of FFTW_R2C_1D_Executor, only the input is n_real_samples/2+1 complex samples.
template <class T>
class FFTW_C2R_1D_Executor {
public:
    typedef fftw_traits<T> traits;

    FFTW_C2R_1D_Executor(int n_real_samples);
    ~FFTW_C2R_1D_Executor();
    void set_input(const typename traits::complex_type* buffer, int size);
    void set_input(const std::vector<typename traits::complex_type>& vec);
    void execute();
    std::vector<T> get_output();

    const int input_size;
    typename traits::complex_type* const input_buffer;

    const int output_size;
    T* const output_buffer;

private:
    typename fftw_traits<T>::plan_type plan;
};

template<class T>
FFTW_R2C_1D_Executor<T>::FFTW_R2C_1D_Executor(int n_real_samples) :
    input_size(n_real_samples),
    input_buffer(fftw_traits<T>::alloc_real(n_real_samples)),
    output_size(n_real_samples/2 + 1),
    output_buffer(fftw_traits<T>::alloc_complex(n_real_samples/2 + 1))
{
    plan = fftw_traits<T>::plan_dft_r2c_1d(n_real_samples, input_buffer, output_buffer, FFTW_ESTIMATE);
}

template<class T>
FFTW_R2C_1D_Executor<T>::~FFTW_R2C_1D_Executor()
{
    fftw_traits<T>::destroy_plan(plan);
    fftw_traits<T>::free(input_buffer);
    fftw_traits<T>::free(output_buffer);
}

template<class T>
void FFTW_R2C_1D_Executor<T>::set_input_zeropadded(const T* buffer, int size)
{
    if (size > input_size) {
        std::cerr << "size: " << size << "input_size: " << input_size << std::endl;
    }
    assert(size <= input_size);
    memcpy(input_buffer, buffer, sizeof(T)*size);
    memset(&input_buffer[size], 0, sizeof(T)*(input_size - size));
}

template<class T>
void FFTW_R2C_1D_Executor<T>::set_input_zeropadded(const std::vector<T>& vec)
{
    set_input_zeropadded(&vec[0], vec.size());
}

template<class T>
void FFTW_R2C_1D_Executor<T>::execute()
{
    fftw_traits<T>::execute(plan);
}

template<class T>
std::vector<typename fftw_traits<T>::complex_type> FFTW_R2C_1D_Executor<T>::get_output()
{
    return std::vector<typename fftw_traits<T>::complex_type>(output_buffer, output_buffer + output_size);
}

template<class T>
FFTW_C2R_1D_Executor<T>::FFTW_C2R_1D_Executor(int n_real_samples) : 
    input_size(n_real_samples/2 + 1),
    input_buffer(fftw_traits<T>::alloc_complex(n_real_samples/2 + 1)),
    output_size(n_real_samples),
    output_buffer(fftw_traits<T>::alloc_real(n_real_samples))
{
    plan = fftw_traits<T>::plan_dft_c2r_1d(n_real_samples, input_buffer, output_buffer, FFTW_ESTIMATE);
}

template<class T>
FFTW_C2R_1D_Executor<T>::~FFTW_C2R_1D_Executor()
{
    fftw_traits<T>::destroy_plan(plan);
    fftw_traits<T>::free(input_buffer);
    fftw_traits<T>::free(output_buffer);
}

template<class T>
void FFTW_C2R_1D_Executor<T>::set_input(const typename fftw_traits<T>::complex_type* buffer, int size)
{
    assert(size == input_size);
    memcpy(input_buffer, buffer, sizeof(typename fftw_traits<T>::complex_type)*size);
    memset(&input_buffer[size], 0, sizeof(typename fftw_traits<T>::complex_type)*(input_size - size));
}

template<class T>
void FFTW_C2R_1D_Executor<T>::set_input(const std::vector<typename fftw_traits<T>::complex_type>& vec)
{
    set_input(&vec[0], vec.size());
}

template<class T>
void FFTW_C2R_1D_Executor<T>::execute()
{
    fftw_traits<T>::execute(plan);
}

template<class T>
std::vector<T> FFTW_C2R_1D_Executor<T>::get_output()
{
    return std::vector<T>(output_buffer, output_buffer + output_size);
}
#endif

