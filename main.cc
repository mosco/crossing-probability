#include <cstdio>
#include <cassert>
extern "C" {
    #include <quadmath.h>
}
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <complex.h>
#include <fftw3.h>

#include "fftw_wrappers.hh"
#include "crossing_probability.hh"

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

extern "C" {
    double incbi(double alpha, double beta, double y);
}

pair<vector<double>, vector<double> > compute_cks_bounds(int n, double cks_statistic_value)
{
    vector<double> lower_bounds;
    vector<double> upper_bounds;
    for (int i = 1; i <= n; ++i) {
        lower_bounds.push_back(incbi(i, n-i+1, cks_statistic_value));
        upper_bounds.push_back(incbi(i, n-i+1, 1-cks_statistic_value));
    }
    return pair<vector<double>, vector<double> >(lower_bounds, upper_bounds);
}

template<class T>
void print_vector(const vector<T>& vec)
{
    for (unsigned int i = 0; i < vec.size(); ++i) {
        cout << vec[i] << ", ";
    }
    cout << endl;
}

template<class T>
void print_array(const T* arr, int n)
{
    for (int i = 0; i < n; ++i) {
        cout << static_cast<double>(arr[i]) << ", ";
    }
    cout << endl;
}

void print_double_array(const double* arr, int n)
{
    for (int i = 0; i < n; ++i) {
        cout << arr[i] << ", ";
    }
    cout << endl;
}

void print_float128_array(const __float128* arr, int length)
{
    for (int i = 0; i < length; ++i) {
        cout << static_cast<double>(arr[i]) << ", ";
    }
    cout << endl;
}

vector<double> convolve(const vector<double>& a, const vector<double>& b)
{
    //cout << "a.size(): " << a.size() << " b.size():" << b.size() << endl;
    assert (a.size() == b.size());
    int n = a.size();
    vector<double> res(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j <= i; ++j) {
            sum += a[j]*b[i-j];
        }
        res[i] = sum;
    }
    return res;
}

void print_usage()
{
    cout << "Usage options:" << endl;
    cout << "    crossing_probability file <boundaries-file>\n";
    cout << "    crossing_probability filefft <boundaries-file>\n";
    cout << "    crossing_probability cks <n-samples> <cks-statistic-value>\n";
    cout << "    crossing_probability cksfft <n-samples> <cks-statistic-value>\n";
    cout << endl;
    cout << "Where <boundaries-file> is a text file that contains a line for every uniform order statistic that specifies the lower and upper boundary for that order statistic, separated by a comma. e.g. \"0.2, 0.7\".\n";
}

int handle_cks_command(const char* n_samples_str, const char* cks_statistic_value_str, bool use_fft)
{
    assert(strlen(n_samples_str) != 0);
    assert(strlen(cks_statistic_value_str) != 0);

    char* endptr = NULL;
    long int n;
    n = strtol(n_samples_str, &endptr, 10);
    if (*endptr != '\0')  {
        cout << "Error parsing n-samples. Got: '" << n_samples_str << "'.\n";
        return -1;
    }

    double cks_value = strtod(cks_statistic_value_str, &endptr);
    if ((cks_value == 0.0) && (endptr == cks_statistic_value_str)) {
        cout << "Error parsing cks-statistic-value. Got: " << cks_statistic_value_str << endl;
        return -1;
    }
    if ((cks_value < 0.0) || (cks_value > 1.0)) {
        cout << "cks-statistic-value must be in the range 0 to 1.\n";
        return -1;
    }

    pair<vector<double>, vector<double> > bounds = compute_cks_bounds(n, cks_value);
    //cout << "Lower bounds:\n";
    //print_vector(bounds.first);
    //cout << "Upper bounds:\n";
    //print_vector(bounds.second);

    double probability = -1;
    if (use_fft) {
        probability = crossing_probability_fft(bounds.first, bounds.second);
    } else {
        probability = crossing_probability(bounds.first, bounds.second);
    }
    cout << probability << endl;

    return 0;
}

int handle_file_command(const char* filename, bool use_fft)
{
    try {
        pair<vector<double>, vector<double> > bounds;
        bounds = read_bounds_file(string(filename));
        verify_bounds_are_valid(bounds.first, bounds.second);
        double probability = -1;
        if (use_fft) {
            probability = crossing_probability_fft(bounds.first, bounds.second);
        } else {
            probability = crossing_probability(bounds.first, bounds.second);
        }
        cout << probability << endl;
        return 0;
    } catch (InputFileReadError e) {
        cout << "Error reading file!\n";
        return -1;
    }
}

void print_complex_array(const double complex* arr, int length)
{
    for (int i = 0; i < length; ++i) {
        cout << creal(arr[i]) << "+" << cimag(arr[i]) << "i, ";
    }
    cout << endl;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        print_usage();
        return -1;
    }

    string command = string(argv[1]);
    if (command == "cks") {
        if (argc != 4) {
            cout << "Wrong number of arguments for 'cks' option.\n";
            print_usage();
            return -1;
        }
        return handle_cks_command(argv[2], argv[3], false);
    } else if (command == "cksfft") {
        if (argc != 4) {
            cout << "wrong number of arguments for 'cksfft64' option.\n";
            print_usage();
            return -1;
        }
        return handle_cks_command(argv[2], argv[3], true);
    } else if (command == "file") {
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
