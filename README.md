crossing-probability
====================

Fast C++ programs for computing the boundary crossing probability of Poisson and Binomial stochastic processes.
The main application of this code is the computation of p-values for continuous goodness-of-fit tests. But it also has potential applications in sequential analysis, change-point detection, queuing theory, diffusion, finance, etc.

More details are available at [my homepage](http://www.wisdom.weizmann.ac.il/~amitmo).

# Building

(the following instructions were tested on Linux and Mac OSX only)


## Prerequisites

* GCC or clang compiler.
* The [FFTW3](http://www.fftw.org/) library. [Installation instructions](http://www.fftw.org/download.html).
* Python (for running the waf build system).


## How to build

Simply run
 ```
 ./waf configure
 ./waf build
 ```

 This will build two programs in the ./build directory:
 * **crossprob** implements algorithms described in [Moscovich-Eiger and Nadler 2015](http://arxiv.org/abs/1503.04363) for computing one-sided and two-sided crossing probabilities for the Poisson and Binomial process.
 * **crossprob_mc** estimates crossing probabilities using Monte-Carlo simulations.
 
You may run the tests with
```./waf test```


## Build errors?

If you installed FFTW3 on your system, the compilation should just work. If FFTW3 is not installed system-wide (e.g. because you do not have root privilieges) then before configuring&building you need to:
* Build FFTW3.
* Set the "FFTW3_INCLUDE_DIR_LOCATION" variable in the "wscript" file to wherever the file "fftw3.h" is at.
* Add the directory containing libfftw3 to the path in the environment variable LD_LIBRARY_PATH.


# Usage

Just run **./crossprob** or **./crossprob_mc**. This will print usage instructions.


# Contact

Feel free to ask any questions: moscovich@gmail.com

Would you like to wrap this code for use in R, Python, Matlab or any other language? That would be wonderful and I'd be happy to help.

Amit Moscovich Eiger.
