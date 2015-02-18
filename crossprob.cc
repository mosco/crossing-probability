#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>

#include <cassert>

#include <fftw3.h>

#include "fftw_wrappers.hh"
#include "two_sided_noncrossing_probability.hh"
#include "string_utils.hh"
#include "read_bounds_file.hh"

using namespace std;

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
    }
}
