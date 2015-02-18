#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <limits>

#include <cassert>

#include <fftw3.h>

#include "fftw_wrappers.hh"
#include "two_sided_noncrossing_probability.hh"
#include "string_utils.hh"

using namespace std;

class InputFileReadError {};

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
    
    if ((upper_bound_steps.front() < 0.0) || (upper_bound_steps.back() > 1.0)) {
        throw runtime_error("Upper bound steps must be in the range 0 to 1.");
    }

    if ((lower_bound_steps.front() < 0.0) || (lower_bound_steps.back() > 1.0)) {
        throw runtime_error("Lower bound steps must be in the range 0 to 1.");
    }

    assert(upper_bound_steps.size() >= lower_bound_steps.size());
    for (int i = 0; i < (int)lower_bound_steps.size(); ++i) {
        if (!(upper_bound_steps[i] <= lower_bound_steps[i])) {
            throw runtime_error("Lower bounary must be lower or equal to the upper boundary in the entire interval [0,1].");
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
    cout << "Lower bound steps: (g(t))\n";
    print_vector(lower_bound_steps);

    getline(f, line);
    vector<double> upper_bound_steps = read_comma_delimited_doubles(line);
    cout << "Upper bound steps: (h(t))\n";
    print_vector(upper_bound_steps);

    verify_bounds_are_valid(lower_bound_steps, upper_bound_steps);
    cout << "Bounds are valid.\n";

    return pair<vector<double>, vector<double> >(lower_bound_steps, upper_bound_steps);
}

void print_usage()
{
    cout << "SYNOPSIS\n";
    cout << "    crossing_probability poisson <n> <boundary-functions-file> [--no-fft]\n";
    cout << "    crossing_probability binomial <n> <boundary-functions-file> [--no-fft]\n";
    cout << endl;
    cout << "DESCRIPTION\n";
    cout << "    crossing_probability poisson <n> <boundary-functions-file> [--no-fft]\n";
    cout << "        Computes the probability that g(t) < xi_n(t) < h(t) for all t in [0,1]\n";
    cout << "        where xi_n(t) is a homogeneous Poisson process of intensity n in the interval [0,1].\n";
    cout << endl;
    cout << "    crossing_probability binomial <n> <boundary-functions-file> [--no-fft]\n";
    cout << "        Computes the probability that g(t) < eta_n(t) < h(t) for all t in [0,1]\n";
    cout << "        where eta_n(t) is the Binomial stochastic process with n samples.\n";
    cout << "        This process is the result of drawing n random variables X_1, ..., X_n\n";
    cout << "        uniformly from the interval [0,1] and constructing the cumulative count function\n";
    cout << "            eta_n(t) = number of X_i < t.\n";
    cout << endl;
    cout << "OPTIONS\n";
    cout << "    <n>\n";
    cout << "        In the Poisson case, this is the intensity of the process (i.e. the expectation of xi_n(1)).\n";
    cout << "        In the Binomial case, this is the number of points drawn from [0,1], hence eta_n(1) = 1.\n";
    cout << endl;
    cout << "    <boundary-functions-file>\n";
    cout << "        This file describes the boundary functions g(t) and h(t).\n";
    cout << "        It must contain exactly 2 lines of monotone-increasing comma-separated numbers between 0 and 1\n";
    cout << "        that are the integer-crossing points of the boundary functions.\n";
    cout << endl;
    cout << "        Line 1: the i-th number in this list is inf{t in [0,1] : g(t) >= i}\n";
    cout << "        Line 2: the i-th number in this list is sup{t in [0,1] : h(t) <= i}\n";
    cout << endl;
    cout << "        Example:\n";
    cout << "            0.3, 0.7, 0.9, 1\n";
    cout << "            0, 0, 0.15, 0.5, 0.8\n";
    cout << endl;
    cout << "    --no-fft\n";
    cout << "        Do not perform convolution using the FFT algorithm.\n";
    cout << "        This is typically faster for small values of n or when g(t) and h(t) are close to each other.\n";
}

int handle_command_line_arguments(int argc, char* argv[])
{
    if (!((argc == 4) || (argc == 5))) {
        print_usage();
        throw runtime_error("Expecting 3 or 4 command line arguments!");
    }

    string command = string(argv[1]);
    long n = string_to_long(argv[2]);
    if (n < 0) {
        print_usage();
        throw runtime_error("n must be non-negative!");
    }
    string filename = string(argv[3]);
    bool use_fft = true;
    if (argc == 5) {
        if (string(argv[4]) == "--no-fft") {
            use_fft = false;
        } else {
            print_usage();
            throw runtime_error("If a 5th command line argument is provided, it must be '--no-fft'");
        }
    }

    pair<vector<double>, vector<double> > bounds = read_bounds_file(filename);

    if (command == "poisson") {
        cout << "Poisson crossing probability: " << 1.0 - poisson_process_noncrossing_probability(n, bounds.first, bounds.second, use_fft) << endl;
    } else if (command == "binomial") {
        cout << "Binomial crossing probability: " << 1.0 - binomial_process_noncrossing_probability(n, bounds.first, bounds.second, use_fft) << endl;
    } else {
        print_usage();
        throw runtime_error("Second command line argument must be 'binomial' or 'poisson'");
    }

    return 0;
}
int main(int argc, char* argv[])
{
    try {
        handle_command_line_arguments(argc, argv);
        return 0;
    } catch (runtime_error& e) {
        cout << "runtime_error exception caught:" << endl;
        cout << e.what() << endl;
        return 1;
    } catch (ifstream::failure& e) {
        cout << "ifstream::failure exception caught:" << endl;
        cout << e.what() << endl;
        return 2;
    } catch (InputFileReadError& e) {
        cout << "Error parsing input file!\n";
        return 3;
    }
}
