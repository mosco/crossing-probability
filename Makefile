FFTW_INCLUDE_LOCATION = /home/amitmo/local/include
C = gcc-4.7
CXX = g++-4.7

#FFTW_INCLUDE_LOCATION = /usr/local/include
#C = gcc-4.9
#CXX = g++-4.9

#CXXFLAGS = -g -pg -Wall -O3 -I$(FFTW_INCLUDE_LOCATION)
#LDFLAGS = -g -pg -lfftw3 -lfftw3q -lquadmath 

#CXXFLAGS = -Wall -O3 -march=native -I$(FFTW_INCLUDE_LOCATION)
#LDFLAGS = -lfftw3 -lfftw3q -lquadmath 
LDFLAGS = -lfftw3

CXXFLAGS = -Wall -g -I$(FFTW_INCLUDE_LOCATION)

OBJECTS = main.o crossing_probability.o 

all: main
	
main: $(OBJECTS) fftw_wrappers.hh crossing_probability.hh
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@ 

crossing_probability.o: crossing_probability.cc fftwconvolver.hh fftw_wrappers.hh

clean:
	rm -rf *.o main

