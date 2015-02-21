#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <cassert>
#include <ctime>
#include <cmath>

#include "string_utils.hh"
#include "read_bounds_file.hh"
#include "tinymt64.h"

using namespace std;

void print_array(const double* arr, int n)
{
    for (int i = 0; i < n; ++i) {
        cout << arr[i] << ", ";
    }
    cout << endl;
}

// Random number generator for uniform samples in the range [0,1]
class UniformRNG {
public:
    UniformRNG(uint64_t seed) {
        cout << "Random seed: " << seed << endl;
        tinymt64_init(&random_state, seed);
    }
    double generate() { return tinymt64_generate_double01(&random_state); }
private:
    tinymt64_t random_state;
};

class ExponentialRNG {
public:
    ExponentialRNG(uint64_t seed, double beta) : rng(seed), beta(beta) {}
    double generate() { return -beta * log(rng.generate()); }
private:
    UniformRNG rng;
    double beta;
};

bool does_integer_step_function_cross(const double* steps, size_t num_steps, const vector<double>& g_steps, const vector<double>& h_steps)
{
    if ((num_steps < g_steps.size()) || (num_steps > h_steps.size())) {
        return true;
    }

    assert(h_steps.size() >= g_steps.size());
    assert(num_steps >= g_steps.size());
    for (size_t i = 0; i < g_steps.size(); ++i) {
        if ((steps[i] > g_steps[i]) || (steps[i] < h_steps[i])) {
            return true;
        }
    }
    assert(h_steps.size() >= num_steps);
    for (size_t i = g_steps.size(); i < num_steps; ++i) {
        if (steps[i] < h_steps[i]) {
            return true;
        }
    }

    return false;
}

double does_random_binomial_process_cross(UniformRNG& rng, vector<double>& tmp_buffer, const vector<double>& g_steps, const vector<double>& h_steps)
{
    for (size_t i = 0; i < tmp_buffer.size(); ++i) {
        tmp_buffer[i] = rng.generate();
    }
    sort(tmp_buffer.begin(), tmp_buffer.end());
    return does_integer_step_function_cross(&tmp_buffer[0], tmp_buffer.size(), g_steps, h_steps);
}

double binomial_process_crossing_probability_montecarlo(long n, const vector<double>& g_steps, const vector<double>& h_steps, long num_simulations)
{
    if ((long)h_steps.size() < n) {
        return 1.0;
    }
    if ((long)g_steps.size() > n) {
        return 1.0;
    }
    UniformRNG rng(time(NULL) + (n<<20));

    vector<double> tmp_buffer(n);
    int count_crossings = 0;
    for (int reps = 0; reps < num_simulations; ++reps) {
        count_crossings += does_random_binomial_process_cross(rng, tmp_buffer, g_steps, h_steps);
    }

    return double(count_crossings) / num_simulations;
}

bool does_random_poisson_process_cross(const vector<double>& g_steps, const vector<double>& h_steps, ExponentialRNG& exprng, vector<double>& tmp_buffer)
{
    size_t max_steps = h_steps.size();
    assert(tmp_buffer.size() >= max_steps);

    size_t num_steps = 0;
    double last_x = 0.0;
    while (true) {
        last_x += exprng.generate();
        if (last_x > 1.0) {
            //cout << "Calling does_integer_step_function_cross with buffer:" << endl;
            //print_array(&tmp_buffer[0], num_steps);
            return does_integer_step_function_cross(&tmp_buffer[0], num_steps, g_steps, h_steps);
        }
        if (num_steps >= max_steps) {
            //cout << "Buffer passed max_steps: (last_x=" << last_x << ")\n";
            //print_array(&tmp_buffer[0], num_steps);
            return true;
        }
        tmp_buffer[num_steps] = last_x;
        ++num_steps;
    }
}

double poisson_process_crossing_probability_montecarlo(double intensity, const vector<double>& g_steps, const vector<double>& h_steps, long num_simulations)
{
    ExponentialRNG exprng(time(NULL) + (int)(intensity*1000000.0), 1.0/intensity);
    vector<double> buffer(h_steps.size() + 1);
    int count_crossings = 0;
    for (int reps = 0; reps < num_simulations; ++reps) {
        count_crossings += does_random_poisson_process_cross(g_steps, h_steps, exprng, buffer);
    }

    return double(count_crossings) / num_simulations;
}

void print_usage()
{
    cout << "SYNOPSIS\n";
    cout << "    crossing_probability poisson <n> <boundary-functions-file> <num-simulations>\n";
    cout << "    crossing_probability binomial <n> <boundary-functions-file> <num-simulations>\n";
    cout << endl;
    cout << "DESCRIPTION\n";
    cout << "    crossing_probability poisson <n> <boundary-functions-file> <num-simulations>\n";
    cout << "        Estimates (using Monte-Carlo simulations) the probability that g(t) < xi_n(t) < h(t) for all t in [0,1]\n";
    cout << "        where xi_n(t) is a homogeneous Poisson process of intensity n in the interval [0,1].\n";
    cout << endl;
    cout << "    crossing_probability binomial <n> <boundary-functions-file> <num-simulations>\n";
    cout << "        Estimates (using Monte-Carlo simulations) the probability that g(t) < eta_n(t) < h(t) for all t in [0,1]\n";
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
    cout << "    <num-simulations>\n";
    cout << "        Number of Monte-Carlo simulation runs.\n";
}

int handle_command_line_arguments(int argc, char* argv[])
{
    if (argc != 5) {
        print_usage();
        throw runtime_error("Expecting 4 command line arguments!");
    }

    string command = string(argv[1]);
    long n = string_to_long(argv[2]);
    if (n < 0) {
        print_usage();
        throw runtime_error("n must be non-negative!");
    }
    string filename = string(argv[3]);
    long num_simulations = string_to_long(argv[4]);
    if (num_simulations < 0) {
        print_usage();
        throw runtime_error("num-simulations must be non-negative!");
    }

    pair<vector<double>, vector<double> > bounds = read_bounds_file(filename);

    if (command == "poisson") {
        cout << "Running " << num_simulations << " simulations...\n";
        double crossprob = poisson_process_crossing_probability_montecarlo(n, bounds.first, bounds.second, num_simulations);
        cout << "Crossing probability: " << crossprob << endl;
    } else if (command == "binomial") {
        cout << "Running " << num_simulations << " simulations...\n";
        double crossprob = binomial_process_crossing_probability_montecarlo(n, bounds.first, bounds.second, num_simulations);
        cout << "Crossing probability: " << crossprob << endl;
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
