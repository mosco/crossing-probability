# Note, FFTW version 3 must be installed to compile this code.
#
# It is assumed that the include file fftw3.h FFTW's static library files are in standard locations that GCC knows about.
# If this is not the case then you need to add the flag -I/wherever/fftw.h/is/located/at to CXXFLAGS
# and add the path of libfftw3.la to the environment variable LIBRARY_PATH.

CXXFLAGS = -Wall -std=c++11 -O3 -march=native -ffast-math
LDFLAGS = -lfftw3

CROSSPROB_OBJECTS = build/crossprob.o build/one_sided_noncrossing_probability.o build/one_sided_noncrossing_probability_n2logn.o build/one_sided_noncrossing_probability_n2.o build/two_sided_noncrossing_probability.o build/fftw_wrappers.o build/fftwconvolver.o build/string_utils.o build/read_boundaries_file.o build/poisson_pmf.o

CROSSPROB_MC_OBJECTS = build/crossprob_mc.o build/string_utils.o build/read_boundaries_file.o build/tinymt64.o


all: crossprob crossprob_mc

crossprob: build bin $(CROSSPROB_OBJECTS)
	$(CXX) $(CROSSPROB_OBJECTS) $(LDFLAGS) -o bin/$@ 

crossprob_mc: build bin $(CROSSPROB_MC_OBJECTS)
	$(CXX) $(CROSSPROB_MC_OBJECTS) $(LDFLAGS) -o bin/$@ 

clean:
	rm -rf build
	rm -rf bin
	rm -rf tests/__pycache__
	rm -rf tests/*.pyc

test: # Running "py.test" also works.
	python tests/test_crossprob.py

swig:
	swig -c++ -python  -o python_extension/crossprob.cc src/crossprob.i

depend:
	makedepend src/*.cc python_extension/*.cc

build/%.o: src/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS)

build:
	mkdir build

bin:
	mkdir bin

# The following dependencies are generated by running "make depend"
# DO NOT DELETE

src/crossprob.o: src/one_sided_noncrossing_probability.hh
src/crossprob.o: src/one_sided_noncrossing_probability_n2logn.hh
src/crossprob.o: src/one_sided_noncrossing_probability_n2.hh
src/crossprob.o: src/two_sided_noncrossing_probability.hh
src/crossprob.o: src/read_boundaries_file.hh src/string_utils.hh
src/crossprob_mc.o: src/string_utils.hh src/read_boundaries_file.hh
src/crossprob_mc.o: src/tinymt64.h /usr/include/stdint.h
src/crossprob_mc.o: /usr/include/features.h /usr/include/stdc-predef.h
src/crossprob_mc.o: /usr/include/inttypes.h
src/fftw_wrappers.o: src/fftw_wrappers.hh /usr/include/fftw3.h
src/fftw_wrappers.o: /usr/include/stdio.h /usr/include/features.h
src/fftw_wrappers.o: /usr/include/stdc-predef.h /usr/include/libio.h
src/fftw_wrappers.o: /usr/include/_G_config.h /usr/include/wchar.h
src/fftw_wrappers.o: src/aligned_mem.hh
src/fftwconvolver.o: src/fftw_wrappers.hh /usr/include/fftw3.h
src/fftwconvolver.o: /usr/include/stdio.h /usr/include/features.h
src/fftwconvolver.o: /usr/include/stdc-predef.h /usr/include/libio.h
src/fftwconvolver.o: /usr/include/_G_config.h /usr/include/wchar.h
src/fftwconvolver.o: src/fftwconvolver.hh src/aligned_mem.hh
src/one_sided_noncrossing_probability.o: src/one_sided_noncrossing_probability.hh
src/one_sided_noncrossing_probability_n2.o: src/one_sided_noncrossing_probability_n2.hh
src/one_sided_noncrossing_probability_n2.o: src/common.hh src/poisson_pmf.hh
src/one_sided_noncrossing_probability_n2.o: src/fftwconvolver.hh
src/one_sided_noncrossing_probability_n2.o: src/aligned_mem.hh
src/one_sided_noncrossing_probability_n2.o: src/string_utils.hh
src/one_sided_noncrossing_probability_n2logn.o: src/one_sided_noncrossing_probability_n2logn.hh
src/one_sided_noncrossing_probability_n2logn.o: src/common.hh
src/one_sided_noncrossing_probability_n2logn.o: src/poisson_pmf.hh
src/one_sided_noncrossing_probability_n2logn.o: src/fftwconvolver.hh
src/one_sided_noncrossing_probability_n2logn.o: src/aligned_mem.hh
src/one_sided_noncrossing_probability_n2logn.o: src/string_utils.hh
src/poisson_pmf.o: src/poisson_pmf.hh src/aligned_mem.hh
src/read_boundaries_file.o: src/read_boundaries_file.hh src/string_utils.hh
src/string_utils.o: src/string_utils.hh
src/tinymt64.o: src/tinymt64.h /usr/include/stdint.h /usr/include/features.h
src/tinymt64.o: /usr/include/stdc-predef.h /usr/include/inttypes.h
src/two_sided_noncrossing_probability.o: src/fftwconvolver.hh
src/two_sided_noncrossing_probability.o: src/aligned_mem.hh src/common.hh
src/two_sided_noncrossing_probability.o: src/poisson_pmf.hh
src/two_sided_noncrossing_probability.o: src/string_utils.hh
python_extension/crossprob.o: /usr/include/string.h /usr/include/features.h
python_extension/crossprob.o: /usr/include/stdc-predef.h
python_extension/crossprob.o: /usr/include/xlocale.h /usr/include/math.h
python_extension/crossprob.o: src/two_sided_noncrossing_probability.hh
python_extension/crossprob.o: src/one_sided_noncrossing_probability.hh
python_extension/crossprob.o: src/one_sided_noncrossing_probability_n2logn.hh
python_extension/crossprob.o: src/one_sided_noncrossing_probability_n2.hh
python_extension/crossprob.o: /usr/include/limits.h
