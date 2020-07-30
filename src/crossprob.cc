#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include "one_sided_noncrossing_probability_n2_old.hh"
#include "one_sided_noncrossing_probability_n2logn.hh"
#include "one_sided_noncrossing_probability_n2.hh"
#include "two_sided_noncrossing_probability.hh"
#include "read_boundaries_file.hh"
#include "string_utils.hh"

using namespace std;

static void print_usage()
{
    cout << "SYNOPSIS\n";
    cout << "    crossprob poisson <n> <boundary-functions-file> [--no-fft]\n";
    cout << "    crossprob ecdf2-n2logn <n> <boundary-functions-file> [--no-fft]\n";
    cout << "    crossprob ecdf1-n2 <n> <one-sided-boundary-functions-file>\n";
    cout << "    crossprob ecdf1-n2logn <n> <one-sided-boundary-functions-file>\n";
    cout << endl;
    cout << "DESCRIPTION\n";
    cout << "    Let g(t), h(t):[0,1] -> R be two functions such that g(t) <= h(t). This program computes\n";
    cout << "    the probability that a Poisson process or an empirical CDF will cross these boundaries.\n";
    cout << "    For more details see https://github.com/mosco/crossing-probability\n";
    cout << endl;
    cout << "    crossprob poisson <n> <boundary-functions-file> [--no-fft]\n";
    cout << "        Computes the probability that a homogeneous Poisson process p(t):[0,1]->{0,1,...} of\n";
    cout << "        intensity n will cross the lower boundary g(t) or the upper boundary h(t) at some point t.\n";
    cout << endl;
    cout << "    crossprob ecdf-n2logn <n> <boundary-functions-file> [--no-fft]\n";
    cout << "        Computes the probability that an empirical CDF F^(t) will cross either the lower\n";
    cout << "        boundary g(t) or the upper boundary h(t) at some point t, where F^(t) is the\n";
    cout << "        empirical CDF of n samples drawn uniformly from the interval [0,1].\n";
    cout << "        i.e. Letting X_1, ..., X_n ~ U[0,1], we have F^(t) = (number of X_i <= t) / n\n";
    cout << "    crossprob ecdf-n2 <n> <boundary-functions-file> [--no-fft]\n";
    cout << "        Faster computation of two-sided CDF crossing\n";
    cout << endl;
    cout << "    crossprob ecdf1-n2-old <n> <one-sided-boundary-functions-file>\n";
    cout << "        Computes the probability that an empirical CDF will cross a single boundary.\n";
    cout << "        This works like the ecdf command above, but using either a lower or upper boundary.\n";
    cout << "        This algorithm is based on repeated integration.\n";
    cout << endl; 
    cout << "    crossprob ecdf1-n2logn <n> <one-sided-boundary-functions-file>\n";
    cout << "        One-sided computation based on O(n^2 logn) two-sided poisson crossing method.\n";
    cout << endl; 
    cout << "    crossprob ecdf1-n2 <n> <one-sided-boundary-functions-file>\n";
    cout << "        Clever O(n^2) algorithm based on the O(n^2 logn) method with larger jumps and corrections. Currently broken.\n";
    cout << endl; 
    cout << "OPTIONS\n";
    cout << "    <n>\n";
    cout << "        In the Poisson case, this is the intensity of the process (i.e. the expectation of p(1)).\n";
    cout << "        In the empirical CDF case, this is the number of points drawn from [0,1].\n";
    cout << endl;
    cout << "    <boundary-functions-file>\n";
    cout << "        This file describes the boundary functions g(t) and h(t).\n";
    cout << "        It must contain 2 lines of monotone-increasing comma-separated numbers between 0 and 1\n";
    cout << endl;
    cout << "        POISSON case - these the integer-crossing points of the boundary functions:\n";
    cout << "        Line 1: the i-th number in this list is inf{t in [0,1] : g(t) >= i}\n";
    cout << "        Line 2: the i-th number in this list is sup{t in [0,1] : h(t) <= i}\n";
    cout << endl;
    cout << "        ECDF case - these are the points at which the boundary functions cross i/n for some integer n.\n";
    cout << "        Line 1: the i-th number in this list is inf{t in [0,1] : g(t) >= i/n}\n";
    cout << "        Line 2: the i-th number in this list is sup{t in [0,1] : h(t) <= i/n}\n";
    cout << endl;
    cout << "        Example:\n";
    cout << "            0.3, 0.7, 0.9, 1\n";
    cout << "            0, 0, 0.15, 0.5, 0.8\n";
    cout << endl;
    cout << "    <one-sided-boundary-functions-file\n";
    cout << "        This file is structured like <boundary-functions-file>, but either the first or second line\n";
    cout << "        must be empty.\n";
    cout << endl;
    cout << "    --no-fft\n";
    cout << "        Revert to algorithm [2] instead of [1].\n";
    cout << endl;
    cout << "NOTES\n";
    cout << "    The two-sided crossing commands poisson/ecdf implement the O(n^2 log n) procedure described in [1] \n";
    cout << "    unless the '--no-fft' flag is used, in which case the O(n^3) procedure of [2] is used instead.\n";
    cout << "    For the ecdf_one_sided command, we use a different algorithm described in [3].\n";
    cout << endl;
    cout << "REFERENCES\n";
    cout << "    [1] A. Moscovich, B. Nadler, Fast calculation of boundary crossing probabilities for\n";
    cout << "        Poisson processes (2016), http://arxiv.org/abs/1503.04363\n";
    cout << "    [2] E. Khmaladze, E. Shinjikashvili, Calculation of noncrossing probabilities for Poisson\n";
    cout << "        processes and its corollaries, Adv. Appl. Probab. 33 (2001) 702-716, http://doi.org/10.1239/aap/\n";
    cout << "    [3] A. Moscovich, B. Nadler, C. Spiegelman, On the exact Berk-Jones statistics and their\n";
    cout << "        p-value calculation (2016), http://arxiv.org/abs/1311.3190\n\n";
}

static int handle_command_line_arguments(int argc, char* argv[])
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

    if (command == "poisson") {
        verify_boundary_is_valid(lower_bound_steps);
        verify_boundary_is_valid(upper_bound_steps);
        if (upper_bound_steps.size() == 0) {
            throw runtime_error("Only a lower boundary is specified. The 'poisson' command currently does not support using just a lower boundary. Sorry about that.");
        }
        vector<double> nocross_probabilities = poisson_process_noncrossing_probability(n, lower_bound_steps, upper_bound_steps, use_fft);
        double total_nocross_probability = 0;
        for (int i = lower_bound_steps.size(); i < upper_bound_steps.size()+1; ++i) {
            total_nocross_probability += nocross_probabilities[i];
        }
        cout << 1.0 - total_nocross_probability << endl;
    } else if (command == "ecdf-n2logn") {
        verify_boundary_is_valid(lower_bound_steps);
        verify_boundary_is_valid(upper_bound_steps);
        cout <<  1.0 - ecdf_noncrossing_probability(n, lower_bound_steps, upper_bound_steps, use_fft) << endl;
    } else if (command == "ecdf-n2") {
        cout << "Not implemented yet.\n";
        //verify_boundary_is_valid(lower_bound_steps);
        //verify_boundary_is_valid(upper_bound_steps);
        //cout <<  1.0 - ecdf_noncrossing_probability_n2(n, lower_bound_steps, upper_bound_steps, use_fft) << endl;
    } else if (command == "ecdf1-n2-old") {
        if (use_fft == false) {
            cout << "Warning: --no-fft flag is superfluous when using the 'ecdf_one_sided' command.\n";
        }
        if ((lower_bound_steps.size() > 0) && (upper_bound_steps.size() > 0)) {
            print_usage();
            throw runtime_error("Expecting EITHER a lower or an upper boundary function when using the 'ecdf_one_sided' command.\n");
        }
        if (upper_bound_steps.size() == 0) {
            verify_boundary_is_valid(lower_bound_steps);
            cout << 1.0 - ecdf_lower_noncrossing_probability_n2_old(n, lower_bound_steps) << endl;
        } else {
            assert(lower_bound_steps.size() == 0);
            verify_boundary_is_valid(upper_bound_steps);
            cout << 1.0 - ecdf_upper_noncrossing_probability_n2_old(n, upper_bound_steps) << endl;
        }
    } else if (command == "ecdf1-n2logn") {
        if (use_fft == false) {
            cout << "Warning: --no-fft flag is superfluous when using the 'ecdf1n2logn' command.\n";
        }
        if ((lower_bound_steps.size() > 0) && (upper_bound_steps.size() > 0)) {
            print_usage();
            throw runtime_error("Expecting EITHER a lower or an upper boundary function when using the 'ecdf1n2logn' command.\n");
        }
        if (upper_bound_steps.size() == 0) {
            verify_boundary_is_valid(lower_bound_steps);
            cout << 1.0 - ecdf_lower_noncrossing_probability_n2logn(n, lower_bound_steps) << endl;
        } else {
            assert(lower_bound_steps.size() == 0);
            verify_boundary_is_valid(upper_bound_steps);
            cout << 1.0 - ecdf_upper_noncrossing_probability_n2logn(n, upper_bound_steps) << endl;
        }
    } else if (command == "ecdf1-n2") {
        if (use_fft == false) {
            cout << "Warning: --no-fft flag is superfluous when using the 'ecdf1n2new' command.\n";
        }
        if ((lower_bound_steps.size() > 0) && (upper_bound_steps.size() > 0)) {
            print_usage();
            throw runtime_error("Expecting EITHER a lower or an upper boundary function when using the 'ecdf1n2new' command.\n");
        }
        if (upper_bound_steps.size() == 0) {
            verify_boundary_is_valid(lower_bound_steps);
            cout << 1.0 - ecdf_lower_noncrossing_probability_n2(n, lower_bound_steps) << endl;
        } else {
            assert(lower_bound_steps.size() == 0);
            verify_boundary_is_valid(upper_bound_steps);
            cout << 1.0 - ecdf_upper_noncrossing_probability_n2(n, upper_bound_steps) << endl;
        }
    } else {
        print_usage();
        throw runtime_error("Second command line argument must be 'poisson', 'ecdf' or 'ecdf_one_sided'");
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
