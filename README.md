crossing-probability
====================

Fast C++ programs (with Python interface) for computing the boundary crossing probability of the empirical CDF.

Currently, the main application of this code is to compute p-values for continuous goodness-of-fit tests such as Kolmogorov-Smirnov and their power for detecting specific distributions. These code may be applied to other problems in the fields of sequential analysis, change-point detection, queuing theory, diffusion, finance, etc.  


# Building the binaries

I've tried to build the code on Linux and Mac OSX only. Building on Windows should be possible (e.g. using MinGW or some other GCC installation), but I haven't tried it. It was tested on the gcc and clang compilers.

Prerequisite: The [FFTW3](http://www.fftw.org/) library. [Installation instructions](http://www.fftw.org/download.html).

Simply run
`make`
This will build two programs in the ./bin directory:

* **bin/crossprob** implements algorithms for computing one-sided and two-sided crossing probabilities.
* **bin/crossprob_mc** estimates crossing probabilities using Monte-Carlo simulations.
 
Then run the tests ```make test```.

# Building the Python extension

Run ```make python``` followed by ```python setup.py install``` to build and install the python module ```crossprob``` into your site-packages directory. This uses the standard distutils system. You should then be able to "import crossprob" in your python code.

## Build errors?

If you installed FFTW3 on your system, the compilation should just work. If FFTW3 is not installed system-wide (e.g. because you do not have root privilieges) then before configuring and building you need to:
* Build FFTW3.
* Add -I<FFTW include dir location> (pointing to wherever "fftw3.h" is located) to CXXFLAGS in the Makefile.
* Add the directory containing libfftw3 to the path in the environment variable LD_LIBRARY_PATH (on linux) or DYLD_LIBRARY_PATH (on OSX).


# Usage

Just run **./bin/crossprob** or **./bin/crossprob_mc**. Usage instructions will be displayed.

