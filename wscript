#!/not/executable/python  # This line is only here so that text editors will use Python syntax highlighting.
import os

APPNAME = 'crossing-probability'
VERSION = '1.1'

PROFILER_SUPPORT = False

def configure(config):
    # If the compiler cannot find the fftw3 include files, place the FFTW3 headers somewhere,
    # say at your home directory under local/include and uncomment the following line:
    # config.env.INCLUDES = ['<your home directory>/local/include/']

    config.load('compiler_c')
    config.load('compiler_cxx')
    config.env.CFLAGS = ['-Wall', '-std=c99', '-ffast-math']
    config.env.CXXFLAGS = ['-O3', '-Wall', '-D__STDC_CONSTANT_MACROS', '-std=c++11', '-ffast-math']
    config.env.LINKFLAGS = ['-g']
    config.env.LIB = ['fftw3']

    if PROFILER_SUPPORT:
        for flags in [config.env.CFLAGS, config.env.CXXFLAGS, config.env.LINKFLAGS]:
            flags.append('-pg')

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

