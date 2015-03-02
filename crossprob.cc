#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>

#include <cassert>

#include <fftw3.h>

#include "fftw_wrappers.hh"
#include "one_sided_noncrossing_probability.hh"
#include "two_sided_noncrossing_probability.hh"
#include "string_utils.hh"
#include "read_boundaries_file.hh"

using namespace std;

void print_usage()
{
    cout << "SYNOPSIS\n";
    cout << "    crossprob poisson2 <n> <boundary-functions-file> [--no-fft]\n";
    cout << "    crossprob binomial2 <n> <boundary-functions-file> [--no-fft]\n";
    cout << "    crossprob binomial1 <n> <one-sided-boundary-functions-file>\n";
    cout << endl;
    cout << "DESCRIPTION\n";
    cout << "    crossprob poisson2 <n> <boundary-functions-file> [--no-fft]\n";
    cout << "        Computes the probability that p(t) will exit the boundaries [g(t), h(t)] at some point t in [0,1]\n";
    cout << "        Where p(t) is a homogeneous Poisson process of intensity n in the interval [0,1].\n";
    cout << endl;
    cout << "    crossprob binomial2 <n> <boundary-functions-file> [--no-fft]\n";
    cout << "        Computes the probability that b(t) will exit the boundaries [g(t), h(t)] at some point t in [0,1]\n";
    cout << "        where b(t) is a binomial process with n samples. This process is the result of drawing n\n";
    cout << "        random variables X_1, ..., X_n uniformly from [0,1] setting b(t) := number of X_i <= t.\n";
    cout << endl;
    cout << "    crossprob binomial1 <n> <one-sided-boundary-functions-file>\n";
    cout << "        Computes the probability that a binomial process with n samples will cross a single boundary.\n";
    cout << "        This works like the binomial2 command above, but using either a lower or upper boundary.\n";
    cout << endl; 
    cout << "OPTIONS\n";
    cout << "    <n>\n";
    cout << "        In the Poisson case, this is the intensity of the process (i.e. the expectation of p(1)).\n";
    cout << "        In the Binomial case, this is the number of points drawn from [0,1], hence b(1) is always n.\n";
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
    cout << "    <one-sided-boundary-functions-file\n";
    cout << "        This file is structured like <boundary-functions-file>, but one of the lines must be empty.\n";
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

    pair<vector<double>, vector<double> > bounds = read_boundaries_file(filename);
    const vector<double>& lower_bound_steps = bounds.first;
    const vector<double>& upper_bound_steps = bounds.second;

    if (command == "poisson2") {
        verify_two_sided_boundaries_are_valid(lower_bound_steps, upper_bound_steps);
        cout << 1.0 - poisson_process_noncrossing_probability(n, lower_bound_steps, upper_bound_steps, use_fft, -1) << endl;
    } else if (command == "binomial2") {
        verify_two_sided_boundaries_are_valid(lower_bound_steps, upper_bound_steps);
        cout <<  1.0 - binomial_process_noncrossing_probability(n, lower_bound_steps, upper_bound_steps, use_fft) << endl;
    } else if (command == "binomial1") {
        if (use_fft == false) {
            cout << "Warning: --no-fft flag is superfluous when using the 'binomial1' command.\n";
        }
        if ((lower_bound_steps.size() > 0) && (upper_bound_steps.size() > 0)) {
            print_usage();
            throw runtime_error("Expecting EITHER a lower or an upper boundary function when using the 'binomial1' command.\n");
        }
        if (upper_bound_steps.size() == 0) {
            verify_one_sided_boundary_is_valid(lower_bound_steps);
            cout << 1.0 - binomial_process_lower_noncrossing_probability(n, lower_bound_steps) << endl;
        } else {
            assert(lower_bound_steps.size() == 0);
            verify_one_sided_boundary_is_valid(upper_bound_steps);
            cout << 1.0 - binomial_process_upper_noncrossing_probability(n, upper_bound_steps) << endl;
        }
    } else {
        print_usage();
        throw runtime_error("Second command line argument must be 'binomial1', 'binomial2' or 'poisson2'");
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
