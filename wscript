#!/not/executable/python  # This line is only here so that text editors will use Python syntax highlighting.
import os
os.environ['CC'] = 'gcc'
os.environ['CXX'] = 'g++'

APPNAME = 'crossing-probability'
VERSION = '1.0'

def configure(config):
    config.load('compiler_c')
    config.load('compiler_cxx')
    config.env.CFLAGS = ['-Wall', '-O3', '-std=c99']
    config.env.CXXFLAGS = ['-DNDEBUG', '-Wall', '-O3', '-D__STDC_CONSTANT_MACROS', '-std=c++98']
    config.env.INCLUDES = ['/home/amitmo/local/include', '/usr/local/include']
    config.env.LIB = ['fftw3', 'quadmath']

def options(opt):
    opt.load('compiler_c')
    opt.load('compiler_cxx')

def dist(ctx):
    ctx.algo = 'zip'
    ctx.files = ctx.path.ant_glob(['*.cc', '*.hh', 'wscript'])
    print ctx.files

def build(bld):
    bld.program(
        source       = 'crossprob.cc one_sided_noncrossing_probability.cc two_sided_noncrossing_probability.cc fftw_wrappers.cc fftwconvolver.cc string_utils.cc read_boundaries_file.cc'.split(),
        target       = 'crossprob',
    )
    bld.objects(source = 'tinymt64.c', name = 'tinymt')
    bld.program(
        source = 'crossprob_mc.cc string_utils.cc read_boundaries_file.cc'.split(),
        target = 'crossprob_mc',
        use = 'tinymt',
    )
