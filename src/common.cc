#include "common.hh"

#include <stdexcept>
#include <sstream>

using namespace std;


static bool is_monotone_increasing(const vector<double>& v)
{
    double prev = -numeric_limits<double>::infinity();
    for (int i = 0; i < (int)v.size(); ++i) {
        if (v[i] < prev) {
            return false;
        }
        prev = v[i];
    }
    return true;
}

void check_boundary_vector(string name, int n, const vector<double>& v)
{
    if ((int)v.size() != n) {
        stringstream ss;
        ss << "Expecting " << n << " input bounds " << name << "_1,...," << name << "_" << n << " but got input of length " << v.size() << ".";
        throw runtime_error(ss.str());
    }

    if (!is_monotone_increasing(v)) {
        stringstream ss;
        ss << name << "_1,...," << name << "_" << n << " must be monotone non-decreasing.";
        throw runtime_error(ss.str());
    }

    if ((v.size() > 0) && ((v.front() < 0.0) || (v.back() > 1.0))) {
        stringstream ss;
        ss << name << "_1,...," << name << "_" << n << " must be in the interval [0,1].";
        throw runtime_error(ss.str());
    }
}

void convolve_same_size(int size, const double* src0, const double* src1, double* dest)
{
    for (int j = 0; j < size; ++j) {
        double convolution_at_j = 0.0;
        for (int k = 0; k <= j; ++k) {
            convolution_at_j += src0[k] * src1[j-k];
        }
        dest[j] = convolution_at_j;
    }
}
