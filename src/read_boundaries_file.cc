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

vector<double> read_and_check_next_line(ifstream& f, string name)
{
    string line;

    getline(f, line);
    vector<double> steps = read_comma_delimited_doubles(line);

    if (!is_monotone_increasing(steps)) {
        throw runtime_error(name + " boundaries are not monotone increasing.");
    }

    if ((steps.size() > 0) && ((steps.front() < 0.0) || (steps.back() > 1.0))) {
        throw runtime_error(name + " boundaries must be in the range 0 to 1.");
    }

    return steps;
}

pair<vector<double>, vector<double> > read_and_check_boundaries_file(string filename)
{
    ifstream f(filename);
    if (!f.is_open()) {
        throw runtime_error("Unable to read input file '" + filename + "'");
    }
    f.exceptions(ifstream::failbit | ifstream::badbit);

    vector<double> b = read_and_check_next_line(f, "b_i");
    vector<double> B = read_and_check_next_line(f, "B_i");

    return pair<vector<double>, vector<double> >(b, B);
}

