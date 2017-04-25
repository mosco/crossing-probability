#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <ctime>
#include <cmath>

#include "string_utils.hh"
#include "read_boundaries_file.hh"
#include "tinymt64.h"

using namespace std;

// Random number generator for uniform samples in the range [0,1]
class RandomNumberGenerator {
public:
    RandomNumberGenerator(uint64_t seed) {
        // cout << "Random seed: " << seed << endl;
        tinymt64_init(&random_state, seed);
    }
    double generate_uniform01()
    {
        return tinymt64_generate_double01(&random_state); 
    }
    double generate_exponential(double beta)
    {
        return -beta * log(generate_uniform01()); 
    }
private:
    tinymt64_t random_state;
};

class ExponentialRNG {
public:
    ExponentialRNG(uint64_t seed, double beta) : rng(seed), beta(beta) {}
    double generate()
    {
        return rng.generate_exponential(beta);
    }
private:
    RandomNumberGenerator rng;
    double beta;
};


static bool does_integer_step_function_cross_lower_boundary(const double* steps, size_t num_steps, const vector<double>& g_steps)
{
    if (num_steps < g_steps.size()) {
        return true;
    } else {
        for (size_t i = 0; i < g_steps.size(); ++i) {
            if (steps[i] > g_steps[i]) {
                return true;
            }
        }
    }
    return false;
}

static bool does_integer_step_function_cross_upper_boundary(const double* steps, size_t num_steps, const vector<double>& h_steps)
{
    if (num_steps > h_steps.size()) {
        return true;
    } else {
        for (size_t i = 0; i < num_steps; ++i) {
            if (steps[i] < h_steps[i]) {
                return true;
            }
        }
        return false;
    }
}

static bool does_integer_step_function_cross(const double* steps, size_t num_steps, const vector<double>& g_steps, const vector<double>& h_steps)
{
    if (h_steps.size() == 0) {
        // Special case, no upper bound specified. Check lower bound crossings only.
        return does_integer_step_function_cross_lower_boundary(steps, num_steps, g_steps);
    } else {
        if ((num_steps < g_steps.size()) || (num_steps > h_steps.size())) {
            return true;
        }
        return does_integer_step_function_cross_lower_boundary(steps, num_steps, g_steps) || does_integer_step_function_cross_upper_boundary(steps, num_steps, h_steps);
    }
}

static double does_random_ecdf_cross(const vector<double>& g_steps, const vector<double>& h_steps, RandomNumberGenerator& rng, vector<double>& tmp_buffer)
{
    double last_x = 0.0;
    for (size_t i = 0; i < tmp_buffer.size(); ++i) {
        last_x += rng.generate_exponential(1);
        tmp_buffer[i] = last_x;
    }
    double normalizing_factor = last_x + rng.generate_exponential(1);
    for (size_t i = 0; i < tmp_buffer.size(); ++i) {
        tmp_buffer[i] /= normalizing_factor;
    }
    return does_integer_step_function_cross(&tmp_buffer[0], tmp_buffer.size(), g_steps, h_steps);
}

static double ecdf_crossing_probability_montecarlo(long n, const vector<double>& g_steps, const vector<double>& h_steps, long num_simulations)
{
    if ((long)g_steps.size() > n) {
        return 1.0;
    }
    if ((h_steps.size() > 0) && ((long)h_steps.size() < n)) {
        return 1.0;
    }

    RandomNumberGenerator rng(time(NULL) + (n<<20));

    vector<double> tmp_buffer(n);
    int count_crossings = 0;
    for (int reps = 0; reps < num_simulations; ++reps) {
        count_crossings += does_random_ecdf_cross(g_steps, h_steps, rng, tmp_buffer);
    }

    return double(count_crossings) / num_simulations;
}

static bool does_random_poisson_process_cross(const vector<double>& g_steps, const vector<double>& h_steps, ExponentialRNG& exprng, vector<double>& tmp_buffer)
{
    size_t max_steps = h_steps.size();
    assert(tmp_buffer.size() >= max_steps);

    size_t num_steps = 0;
    double last_x = 0.0;
    while (true) {
        last_x += exprng.generate();
        if (last_x > 1.0) {
            return does_integer_step_function_cross(&tmp_buffer[0], num_steps, g_steps, h_steps);
        }
        if (num_steps >= max_steps) {
            return true;
        }
        tmp_buffer[num_steps] = last_x;
        ++num_steps;
    }
}

static double poisson_process_crossing_probability_montecarlo(double intensity, const vector<double>& g_steps, const vector<double>& h_steps, long num_simulations)
{
    if ((h_steps.size() > 0) && (h_steps.size() < g_steps.size())) {
        return 1.0;
    }

    ExponentialRNG exprng(time(NULL) + (int)(intensity*1000000.0), 1.0/intensity);

    vector<double> buffer(h_steps.size() + 1);
    int count_crossings = 0;
    for (int reps = 0; reps < num_simulations; ++reps) {
        count_crossings += does_random_poisson_process_cross(g_steps, h_steps, exprng, buffer);
    }

    return double(count_crossings) / num_simulations;
}

static void print_usage()
{
    cout << "SYNOPSIS\n";
    cout << "    crossing_probability poisson <n> <boundary-functions-file> <num-simulations>\n";
    cout << "    crossing_probability ecdf <n> <boundary-functions-file> <num-simulations>\n";
    cout << endl;
    cout << "DESCRIPTION\n";
    cout << "    crossing_probability poisson <n> <boundary-functions-file> <num-simulations>\n";
    cout << "        Estimates (using Monte-Carlo simulations) the probability that g(t) < xi_n(t) < h(t) for all t in [0,1]\n";
    cout << "        where xi_n(t) is a homogeneous Poisson process of intensity n in the interval [0,1].\n";
    cout << endl;
    cout << "    crossing_probability ecdf <n> <boundary-functions-file> <num-simulations>\n";
    cout << "        Estimates (using Monte-Carlo simulations) the probability that g(t) < F_n(t) < h(t) for all t in [0,1]\n";
    cout << "        where F_n(t) is the empirical CDF of n uniform samples in [0,1]. i.e.\n";
    cout << "            F_n(t) = (number of X_i < t)/n  where X_1,...X_n ~ U[0,1].\n";
    cout << endl;
    cout << "OPTIONS\n";
    cout << "    <n>\n";
    cout << "        In the Poisson case, this is the intensity of the process (i.e. the expectation of xi_n(1)).\n";
    cout << "        In the empirical CDF case, this is the number of points drawn from [0,1].\n";
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

static int handle_command_line_arguments(int argc, char* argv[])
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

    pair<vector<double>, vector<double> > bounds = read_boundaries_file(filename);
    verify_boundary_is_valid(bounds.first);
    verify_boundary_is_valid(bounds.second);

    if (command == "poisson") {
        // cout << "Running " << num_simulations << " simulations...\n";
        double crossprob = poisson_process_crossing_probability_montecarlo(n, bounds.first, bounds.second, num_simulations);
        cout << crossprob << endl;
    } else if (command == "ecdf") {
        // cout << "Running " << num_simulations << " simulations...\n";
        double crossprob = ecdf_crossing_probability_montecarlo(n, bounds.first, bounds.second, num_simulations);
        cout << crossprob << endl;
    } else {
        print_usage();
        throw runtime_error("Second command line argument must be 'ecdf' or 'poisson'");
    }

    return 0;
}

int main(int argc, char* argv[])
{
    try {
        handle_command_line_arguments(argc, argv);
        return 0;
    } catch (ifstream::failure& e) {
        cout << "ifstream::failure exception caught:" << endl;
        cout << e.what() << endl;
        return 1;
    } catch (runtime_error& e) {
        cout << "runtime_error exception caught:" << endl;
        cout << e.what() << endl;
        return 2;
    }
}
