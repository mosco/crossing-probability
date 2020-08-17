// This is a SWIG interface file for wrapping the C++ code to make it available in other languagses.
//
// The python extension is under python_extension/crossprob_swig_wrap.cc, to rebuild it install SWIG and run "make python".

%define DOCSTRING
"crossprob - crossing probabilities for one-dimensional empirical processes

Let X_1,...,X_n be a set of points sampled uniformly from the interval [0,1]
and let X_(1) <= X_(2) <= ... <= X_(n) be the sorted sample.

This module implements several algorithms for computing the probability
    Pr[for all i: b_i <= X_(i) <= B_i]        (1)

Equivalently, let the empirical cumulative distribution function of X_1,...,X_n be
    F_n(t) := \sum_i 1(X_i <= t).
The probability in Eq. (1) is equivalent to the probability that F_n(t) does not cross
a two-sided boundary,

The most versatile function is the two-sided crossing probability
    ecdf2(b, B, use_fft)
Its parameters are:
    n: the size of the sample
    b, B: two lists of length n of the boundaries in Eq. (1) above.
    use_fft: If true algorithm the O(n^2 logn) algorithm [MNS2016] is used,
             otherwise the O(n^3) algorithm of [KS2001]

Faster functions are available for the special case of a single boundary:
    ecdf1_new_b(b)
        Implements a new O(n^2) algorithm. B_i are implicitly assumed to be 1. 
    ecdf1_new_B(B)
        Implements the O(n^2) algorithm. b_i are implicitly assumed to be 0. 
    ecdf1_mns2016_b(b)
        Implements the O(n^2) algorithm of [MNS2016]. B_i are implicitly assumed to be 1. 
        Generally slower and less numerically stable than ecdf1_new_b()
    ecdf1_mns2016_B(B)
        Implements the O(n^2) algorithm of [MNS2016]. b_i are implicitly assumed to be 0. 
        Generally slower and less numerically stable than ecdf1_new_B()

EXAMPLES
    For a sample X_1, X_2, X_3 with order statistics X_(1) <= X_(2) <= X(3), the probability
        Pr[X_(1)<=0.7 and 0.15<=X_(2)<=0.9 and 0.5<=X_(3)<=0.8]
    may be computed using
        ecdf2([0, 0.15, 0.5], [0.7, 0.9, 0.8], True)

REFERENCES
    [KS2001] Estate Khmaladze, Eka Shinjikashvili (2001). Calculation of noncrossing probabilities for Poisson
             processes and its corollaries, Advances in Applied Probability. https://doi.org/10.1239/aap/1005091361
    [MNS2016] Amit Moscovich, Boaz Nadler, Clifford Spiegelman (2016). On the exact Berk-Jones statistics and their
              p-value calculation. Electronic Journal of Statistics. https://doi.org/10.1214/16-EJS1172
    [MN2017] Amit Moscovich, Boaz Nadler (2017). Fast calculation of boundary crossing probabilities for Poisson processes.
             Statistics & Probability Letters. https://doi.org/10.1016/j.spl.2016.11.027

MODULE REFERENCE
    https://github.com/mosco/crossing-probability
"
%enddef

%module(docstring=DOCSTRING) crossprob


%include "std_vector.i"
namespace std {
   %template(VectorDouble) vector<double>;
};

%exception {
    try {
        $action
    } catch(std::runtime_error& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch(std::exception& e) {
        SWIG_exception(SWIG_UnknownError, e.what());
    } catch(...) {
        SWIG_exception(SWIG_UnknownError, "Caught unknown C++ exception (not an std::exception object)");
    }
}

%{
#include "../src/ecdf2.hh"
#include "../src/ecdf1_mns2016.hh"
#include "../src/ecdf1_new.hh"
%}

%feature("autodoc", "1");
%include "../src/ecdf2.hh"
%include "../src/ecdf1_mns2016.hh"
%include "../src/ecdf1_new.hh"

