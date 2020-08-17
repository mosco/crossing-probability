// This is a SWIG interface file for wrapping the C++ code to make it available in other languagses.
//
// The python extension is under python_extension/crossprob_swig_wrap.cc, to rebuild it install SWIG and run "make python".

%module crossprob

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
#include "../src/one_sided_noncrossing_probability_n2_old.hh"
#include "../src/one_sided_noncrossing_probability_n2.hh"
#include "../src/two_sided_noncrossing_probability.hh"
%}

%include "../src/one_sided_noncrossing_probability_n2_old.hh"
%include "../src/one_sided_noncrossing_probability_n2.hh"
%include "../src/two_sided_noncrossing_probability.hh"

