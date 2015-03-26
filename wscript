#!/not/executable/python  # This line is only here so that text editors will use Python syntax highlighting.
import os


# If the compiler cannot find the fftw3 include files, add a line like the following line to the configure() function:
#     config.env.INCLUDES = ['/location/of/fftw3/include/directory']
#
# To override the selected compiler, uncomment the following lines and choose a compiler:
#     os.environ['CC'] = 'gcc-4.9'
#     os.environ['CXX'] = 'g++-4.9'

APPNAME = 'crossing-probability'
VERSION = '1.0'

def configure(config):
    config.load('compiler_c')
    config.load('compiler_cxx')
    config.env.CFLAGS = ['-Wall', '-O3', '-std=c99']
    config.env.CXXFLAGS = ['-DNDEBUG', '-Wall', '-O3', '-D__STDC_CONSTANT_MACROS', '-std=c++11']
    config.env.LIB = ['fftw3']

def options(opt):
    opt.load('compiler_c')
    opt.load('compiler_cxx')

def dist(ctx):
    ctx.algo = 'zip'
    ctx.files = ctx.path.ant_glob(['src/*', 'tests/*.txt', 'tests/*.py', 'wscript', 'waf', 'README.md'])
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

