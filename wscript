#!/not/executable/python  # This line is only here so that text editors will use Python syntax highlighting.
import os

C_COMPILER = 'gcc-4.7'
CPLUSPLUS_COMPILER = 'g++-4.7'
FFTW3_INCLUDE_DIR_LOCATION = '/home/amitmo/local/include' 
# Note: the libfftw3 file should be in a path contained in LD_LIBRARY_PATH.

os.environ['CC'] = C_COMPILER
os.environ['CXX'] = CPLUSPLUS_COMPILER

APPNAME = 'crossing-probability'
VERSION = '1.0'

def configure(config):
    config.load('compiler_c')
    config.load('compiler_cxx')
    config.env.CFLAGS = ['-Wall', '-O3', '-std=c99']
    config.env.CXXFLAGS = ['-DNDEBUG', '-Wall', '-O3', '-D__STDC_CONSTANT_MACROS', '-std=c++98']
    config.env.INCLUDES = [FFTW3_INCLUDE_DIR_LOCATION]
    config.env.LIB = ['fftw3']

def options(opt):
    opt.load('compiler_c')
    opt.load('compiler_cxx')

def dist(ctx):
    ctx.algo = 'zip'
    ctx.files = ctx.path.ant_glob(['src/*', 'wscript', 'waf', 'README.md'])
    print ctx.files

def make_cc_file_list(cc_files_string):
    names = cc_files_string.split()
    return [('src/%s.cc' % name) for name in names]

def build(bld):
    bld.program(
        source       = make_cc_file_list('''
            crossprob
            one_sided_noncrossing_probability
            two_sided_noncrossing_probability
            fftw_wrappers
            fftwconvolver
            string_utils
            read_boundaries_file
        '''),
        target       = 'crossprob',
    )

    bld.objects(source = 'src/tinymt64.c', name = 'tinymt')

    bld.program(
        use = 'tinymt',
        source = make_cc_file_list('''
            crossprob_mc
            string_utils
            read_boundaries_file
        '''),
        target = 'crossprob_mc',
    )

def test(ctx):
    import tests.test_crossprob
    tests.test_crossprob.main()

