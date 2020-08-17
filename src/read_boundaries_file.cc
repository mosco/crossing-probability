#include "read_boundaries_file.hh"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cassert>
#include <limits>

#include "string_utils.hh"

using namespace std;


vector<double> read_next_line(ifstream& f, string name)
{
    string line;

    getline(f, line);
    vector<double> steps = read_comma_delimited_doubles(line);

    return steps;
}

pair<vector<double>, vector<double> > read_and_check_boundaries_file(string filename)
{
    ifstream f(filename);
    if (!f.is_open()) {
        throw runtime_error("Unable to read input file '" + filename + "'");
    }
    f.exceptions(ifstream::failbit | ifstream::badbit);

    vector<double> b = read_next_line(f, "b_i");
    vector<double> B = read_next_line(f, "B_i");

    return pair<vector<double>, vector<double> >(b, B);
}

