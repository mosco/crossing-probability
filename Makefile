FFTW_INCLUDE_LOCATION = /home/amitmo/local/include
C = gcc-4.7
CXX = g++-4.7

#FFTW_INCLUDE_LOCATION = /usr/local/include
#C = gcc-4.9
#CXX = g++-4.9

#CXXFLAGS = -g -pg -Wall -O3 -I$(FFTW_INCLUDE_LOCATION)
#LDFLAGS = -g -pg -lfftw3 -lfftw3q -lquadmath 

CXXFLAGS = -Wall -O3 -march=native -I$(FFTW_INCLUDE_LOCATION)
LDFLAGS = -lfftw3

all: crossprob 
	
OBJECTS_CROSSPROB = crossprob.o twosided_noncrossing_probability.o fftw_wrappers.o fftwconvolver.o string_utils.o
crossprob: $(OBJECTS_CROSSPROB) fftw_wrappers.hh twosided_noncrossing_probability.hh string_utils.hh
	$(CXX) $(OBJECTS_CROSSPROB) $(LDFLAGS) -o $@ 

twosided_noncrossing_probability.o: twosided_noncrossing_probability.hh fftwconvolver.hh fftw_wrappers.hh
fftw_wrappers.o: fftw_wrappers.hh
fftwconvolver.o: fftw_wrappers.hh fftwconvolver.hh

clean:
	rm -rf *.o crossprob

