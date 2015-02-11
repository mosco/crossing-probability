#include <cassert>
#include <iostream>
#include <vector>
#include <fftw3.h>

#include "fftw_wrappers.hh"
#include "noncrossing_probability.hh"

using namespace std;

class InputFileReadError {};

pair<vector<double>, vector<double> > read_bounds_file(const string filename)
{
    vector<double> lower_bounds, upper_bounds;
    FILE* fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        throw new InputFileReadError();
    }

    int line_number = 1;
    while (true) {
        double lower_bound;
        double upper_bound;
        int ret;
        ret = fscanf(fp, "%lf , %lf", &lower_bound, &upper_bound);
        if (ret == 2) {
            lower_bounds.push_back(lower_bound);
            upper_bounds.push_back(upper_bound);
        } else if (ret == EOF) {
            break;
        } else {
            cout << "Error on line " << line_number << ":" << endl;
            cout << "Unable to parse line.";
            throw InputFileReadError();
        }
        ++line_number;
    }

    return pair<vector<double>, vector<double> >(lower_bounds, upper_bounds);
}

void verify_bounds_are_valid(const vector<double>& lower_bounds, const vector<double>& upper_bounds)
{
    assert (lower_bounds.size() == upper_bounds.size());
    int n = lower_bounds.size();
    double lower_bound_min_value = 0.0;
    double upper_bound_min_value = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double lower_bound = lower_bounds[i];
        if ((lower_bound < lower_bound_min_value) || (lower_bound > 1.0)) {
            cout << "Error on line " << i+1 << ":" << endl;
            cout << "Lower bounds must form an increasing sequence of numbers between 0 and 1.\n";
            cout << "Expected lower bound L to satisfy " << lower_bound_min_value << " <= L <= 1.0. Instead got " << lower_bound << ".\n";
            throw InputFileReadError();
        }
        lower_bound_min_value = lower_bound;

        double upper_bound = upper_bounds[i];
        if ((upper_bound < upper_bound_min_value) || (upper_bound > 1.0)) {
            cout << "Error on line " << i+1 << ":" << endl;
            cout << "Lower bounds must form an increasing sequence of numbers between 0 and 1.\n";
            cout << "Expected upper bound U to satisfy " << lower_bound_min_value << " <= U <= 1.0. Got " << upper_bound << ".\n";
            throw InputFileReadError();
        }

        if (lower_bound > upper_bound) {
            cout << "Error on line " << i+1 << ":" << endl;
            cout << "Lower bound is larger than upper bound. Lower bound: " << lower_bound << ", upper bound: " << upper_bound << endl;
        }
        upper_bound_min_value = upper_bound;
    }
}

void print_usage()
{
    cout << "Usage options:" << endl;
    cout << "    crossing_probability file <boundaries-file>\n";
    cout << "    crossing_probability filefft <boundaries-file>\n";
    cout << endl;
    cout << "Where <boundaries-file> is a text file that contains a line for every uniform order statistic that specifies the lower and upper boundary for that order statistic, separated by a comma. e.g. \"0.2, 0.7\".\n";
}

int handle_file_command(const char* filename, bool use_fft)
{
    try {
        pair<vector<double>, vector<double> > bounds;
        bounds = read_bounds_file(string(filename));
        verify_bounds_are_valid(bounds.first, bounds.second);
        double probability = -1;
        if (use_fft) {
            probability = 1.0 - binomial_process_noncrossing_probability_fft(bounds.first, bounds.second);
        } else {
            probability = 1.0 - binomial_process_noncrossing_probability(bounds.first, bounds.second);
        }
        cout << probability << endl;
        return 0;
    } catch (InputFileReadError e) {
        cout << "Error reading file!\n";
        return -1;
    }
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        print_usage();
        return -1;
    }

    string command = string(argv[1]);
    if (command == "file") {
        if (argc != 3) {
            cout << "Wrong number of arguments for 'file' option.\n";
            print_usage();
            return -1;
        }
        return handle_file_command(argv[2], false);
    } else if (command == "filefft") {
        if (argc != 3) {
            cout << "Wrong number of arguments for 'file' option.\n";
            print_usage();
            return -1;
        }
        return handle_file_command(argv[2], true);
    } else {
        cout << "Unknown option '" << argv[1] << "'.\n";
        return -1;
    }
}
