#include "read_bounds_file.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <stdexcept>
#include <cassert>

#include "string_utils.hh"

using namespace std;

bool is_monotone_increasing(const vector<double>& v)
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

void verify_bounds_are_valid(const vector<double>& lower_bound_steps, const vector<double>& upper_bound_steps)
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

template <class T>
void print_vector(const vector<T>& v)
{
    for (int i = 0; i < (int)v.size(); ++i) {
        cout << v[i] << ", ";
    }
    cout << endl;
}

pair<vector<double>, vector<double> > read_bounds_file(const string& filename)
{
    string line;
    ifstream f(filename.c_str());
    if (!f.is_open()) {
        throw runtime_error("Unable to read input file '" + filename + "'");
    }
    f.exceptions(ifstream::failbit | ifstream::badbit);

    getline(f, line);
    vector<double> lower_bound_steps = read_comma_delimited_doubles(line);
    //cout << "Lower bound steps: (g(t))\n";
    //print_vector(lower_bound_steps);

    getline(f, line);
    vector<double> upper_bound_steps = read_comma_delimited_doubles(line);
    //cout << "Upper bound steps: (h(t))\n";
    //print_vector(upper_bound_steps);

    verify_bounds_are_valid(lower_bound_steps, upper_bound_steps);
    cout << "Bounds are valid.\n";

    return pair<vector<double>, vector<double> >(lower_bound_steps, upper_bound_steps);
}

