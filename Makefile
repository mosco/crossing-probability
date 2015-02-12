#FFTW_INCLUDE_LOCATION = /home/amitmo/local/include
#C = gcc-4.7
#CXX = g++-4.7

FFTW_INCLUDE_LOCATION = /usr/local/include
C = gcc-4.9
CXX = g++-4.9

#CXXFLAGS = -g -pg -Wall -O3 -I$(FFTW_INCLUDE_LOCATION)
#LDFLAGS = -g -pg -lfftw3 -lfftw3q -lquadmath 

CXXFLAGS = -Wall -O3 -march=native -I$(FFTW_INCLUDE_LOCATION)
#LDFLAGS = -lfftw3 -lfftw3q -lquadmath 
LDFLAGS = -lfftw3

#CXXFLAGS = -Wall -g -I$(FFTW_INCLUDE_LOCATION)

OBJECTS = main.o noncrossing_probability.o fftw_wrappers.o fftwconvolver.o

all: main
	
main: $(OBJECTS) fftw_wrappers.hh noncrossing_probability.hh
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@ 

noncrossing_probability.o: noncrossing_probability.cc noncrossing_probability.hh fftwconvolver.hh fftw_wrappers.hh
fftw_wrappers.o: fftw_wrappers.hh
fftwconvolver.o: fftw_wrappers.hh fftwconvolver.hh

clean:
	rm -rf *.o main

