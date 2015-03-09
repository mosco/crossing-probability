#include "read_boundaries_file.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <limits>

#include "string_utils.hh"

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

pair<vector<double>, vector<double> > read_boundaries_file(const string& filename)
{
    string line;
    ifstream f(filename.c_str());
    if (!f.is_open()) {
        throw runtime_error("Unable to read input file '" + filename + "'");
    }
    f.exceptions(ifstream::failbit | ifstream::badbit);

    getline(f, line);
    vector<double> lower_bound_steps = read_comma_delimited_doubles(line);

    getline(f, line);
    vector<double> upper_bound_steps = read_comma_delimited_doubles(line);

    return pair<vector<double>, vector<double> >(lower_bound_steps, upper_bound_steps);
}

void verify_one_sided_boundary_is_valid(const vector<double>& steps)
{
    if (!is_monotone_increasing(steps)) {
        throw runtime_error("Bound steps are not monotone increasing.");
    }
    if ((steps.size() > 0) && ((steps.front() < 0.0) || (steps.back() > 1.0))) {
        throw runtime_error("Steps must be in the range 0 to 1.");
    }
}

void verify_two_sided_boundaries_are_valid(const vector<double>& lower_bound_steps, const vector<double>& upper_bound_steps)
{
    int n = upper_bound_steps.size();
    if ((int)lower_bound_steps.size() > n) {
        throw runtime_error("The number of lower bound steps must be smaller or equal to the number of upper bound steps, otherwise this means the lower boundary crosses the upper boundary.");
    }

    if (!is_monotone_increasing(upper_bound_steps)) {
        throw runtime_error("Upper bound steps are not monotone increasing.");
    }

    if (!is_monotone_increasing(lower_bound_steps)) {
        throw runtime_error("Lower bound steps are not monotone increasing."); 
    }
    
    if ((upper_bound_steps.size() > 0) && ((upper_bound_steps.front() < 0.0) || (upper_bound_steps.back() > 1.0))) {
        throw runtime_error("Upper bound steps must be in the range 0 to 1.");
    }

    if ((lower_bound_steps.size() > 0) && ((lower_bound_steps.front() < 0.0) || (lower_bound_steps.back() > 1.0))) {
        throw runtime_error("Lower bound steps must be in the range 0 to 1.");
    }

    assert(upper_bound_steps.size() >= lower_bound_steps.size());
    for (int i = 0; i < (int)lower_bound_steps.size(); ++i) {
        if (!(upper_bound_steps[i] <= lower_bound_steps[i])) {
            stringstream ss;
            ss << "Lower boundary must be lower or equal to the upper boundary in the entire interval [0,1]: upper_bound_steps[" << i << "]=" << upper_bound_steps[i] << " lower_bound_steps[" << i << "]=" << lower_bound_steps[i];
            throw runtime_error(ss.str());
        }
    }
}

