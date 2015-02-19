#FFTW_INCLUDE_LOCATION = /home/amitmo/local/include
#C = gcc-4.7
#CXX = g++-4.7

FFTW_INCLUDE_LOCATION = /usr/local/include
C = gcc-4.9
CXX = g++-4.9

#CXXFLAGS = -g -pg -Wall -O3 -I$(FFTW_INCLUDE_LOCATION)
#LDFLAGS = -g -pg -lfftw3 -lfftw3q -lquadmath 

CXXFLAGS = -Wall -O3 -march=native -I$(FFTW_INCLUDE_LOCATION)
LDFLAGS = -lfftw3

TARGETS = crossprob crossprobmc

all: $(TARGETS)
	
OBJECTS_CROSSPROB = crossprob.o two_sided_noncrossing_probability.o fftw_wrappers.o fftwconvolver.o string_utils.o read_bounds_file.o
crossprob: $(OBJECTS_CROSSPROB) fftw_wrappers.hh two_sided_noncrossing_probability.hh string_utils.hh read_bounds_file.hh
	$(CXX) $(OBJECTS_CROSSPROB) $(LDFLAGS) -o $@ 

OBJECTS_CROSSPROB_MC = crossprob_mc.o string_utils.o read_bounds_file.o
crossprobmc: $(OBJECTS_CROSSPROB_MC) string_utils.hh read_bounds_file.hh
	$(CXX) $(OBJECTS_CROSSPROB_MC) $(LDFLAGS) -o $@  

%.o: %.cc Makefile
	@$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

#two_sided_noncrossing_probability.o: two_sided_noncrossing_probability.hh fftwconvolver.hh fftw_wrappers.hh
#fftw_wrappers.o: fftw_wrappers.hh
#fftwconvolver.o: fftw_wrappers.hh fftwconvolver.hh
#string_utils.o: string_utils.hh
#read_bounds_file.o: read_bounds_file.hh string_utils.hh

clean:
	rm -rf *.o $(TARGETS)

