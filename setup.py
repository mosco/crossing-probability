from distutils.core import setup, Extension
import os

module1 = Extension(
    '_crossprob',
    sources = [
        'src/one_sided_noncrossing_probability_n2_old.cc',
        'src/one_sided_noncrossing_probability_n2.cc',
        'src/two_sided_noncrossing_probability.cc',
        'src/fftw_wrappers.cc',
        'src/fftwconvolver.cc',
        'src/string_utils.cc',
        'src/poisson_pmf.cc',
        'python_extension/crossprob.cc'
    ],
    extra_compile_args = ['-Wall', '-std=c++11', '-O3', '-march=native', '-ffast-math'],
    extra_link_args = ['-lfftw3'],
    undef_macros = ["NDEBUG"]
)

setup(name = 'crossprob',
    version = '1.2',
    description = 'Fast computation of boundary crossing probabilities for Poisson and empirical processes, for more details see https://github.com/mosco/crossing-probability',
    ext_modules = [module1],
    py_modules = ['crossprob'],
    package_dir={'': 'python_extension'}
)
