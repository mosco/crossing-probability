// This is a SWIG interface file for wrapping the C++ code to make it available in other languagses.
//
// The python extension is under python_extension/crossprob_swig_wrap.cc, to rebuild it install SWIG and run "make swig".

%module crossprob

%include "std_vector.i"
namespace std {
   %template(VectorDouble) vector<double>;
};

%{
#include "../src/two_sided_noncrossing_probability.hh"
#include "../src/one_sided_noncrossing_probability.hh"
#include "../src/one_sided_noncrossing_probability_n2logn.hh"
%}

%include "../src/one_sided_noncrossing_probability.hh"
%include "../src/two_sided_noncrossing_probability.hh"
%include "../src/one_sided_noncrossing_probability_n2logn.hh"


