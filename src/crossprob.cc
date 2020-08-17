#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <algorithm>

#include "read_boundaries_file.hh"
#include "string_utils.hh"

#include "ecdf1_mns2016.hh"
#include "ecdf1_new.hh"
#include "ecdf2.hh"

using namespace std;

static void print_usage()
{
    cout << "SYNOPSIS\n";
    cout << "    crossprob <algorithm> <n> <one-sided-boundary-file>\n";
    cout << "\n";
    cout << "DESCRIPTION\n";
    cout << "    Let X_1, ..., X_n be a set of points sampled uniformly from the interval [0,1]\n";
    cout << "    and let X_(1) <= X_(2) <= ... <= X_(n) be the sorted sample.\n";
    cout << "\n";
    cout << "    This program implements several algorithms for computing the probability that\n";
    cout << "        for all i: b_i <= X_(i) <= B_i\n"; 
    cout << "    It also has one-sided crossing variants:\n";
    cout << "        for all i: b_i <= X_(i)\n";
    cout << "    and\n";
    cout << "        for all i: X_(i) <= B_i\n";
    cout << "\n";
    cout << "    For more details see https://github.com/mosco/crossing-probability\n";
    cout << "    and the references at the end of this help page.\n";
    cout << "\n";
    cout << "OPTIONS\n";
    cout << "    <algorithm>\n";
    cout << "        ecdf2-ks2001: an O(n^3) algorithm for two-sided boundaries. [1]\n";
    cout << "        ecdf2-mn2017: an O(n^2 log n) method for two-sided boundaries. [3]\n";
    cout << "        ecdf1-mns2016: an O(n^2) method for one-sided boundaries. [2]\n";
    cout << "        ecdf1-new: New O(n^2) method, typically faster than ecdf1-mns2016.\n";
    cout << "\n";            
    cout << "    <n>\n";
    cout << "        The number of sample points drawn from [0,1].\n";
    cout << "\n";
    cout << "    <boundary-functions-file>\n";
    cout << "        This text file contains the two lines of comma-separater numbers:\n";
    cout << "            b_1, b_2, ..., b_n\n";
    cout << "            B_1, B_2, ..., B_n\n";
    cout << "\n";
    cout << "        For the two-sided ecdf2-* algorithms, both lines should have n numbers.\n";
    cout << "        If the first line is empty b_i is implicitly assumed to be 0.\n";
    cout << "        If the second line is empty B_i is implicitly assumed to be 1.\n";
    cout << "\n";
    cout << "        For the one-sided ecdf1-* algorithms one of the input lines must\n";
    cout << "        have n elements and the other line must be empty.\n";
    cout << "\n";
    cout << "EXAMPLES:\n";
    cout << "    To check the probability that\n";
    cout << "    X_(1)<=0.7 and 0.15<=X_(2)<=0.9 and 0.5<=X_(3)<= 0.7\n";
    cout << "    write the following bounds.txt file:\n";
    cout << "    0, 0.15, 0.5,\n";
    cout << "    0.7, 0.9, 1\n";
    cout << "    and run './bin/crossprob ecdf2-mn2017 3 bounds.txt'\n";
    cout << "\n";
    cout << "    To compute a one-sided crossing probability for two samples\n";
    cout << "    that X_(1) <= 0.5 and X_(2) <= 0.7, we can run\n";
    cout << "        ./bin/crossprob ecdf1-m2020 2 bounds1.txt\n";
    cout << "    where bounds1.txt is the following (first line is empty):\n";
    cout << "            \n";
    cout << "    0.5, 0.7\n";
    cout << "\n";
    cout << "REFERENCES\n";
    cout << "    [1] Estate Khmaladze, Eka Shinjikashvili (2001). Calculation of noncrossing probabilities for Poisson\n";
    cout << "        processes and its corollaries, Advances in Applied Probability. https://doi.org/10.1239/aap/1005091361\n";
    cout << "    [2] Amit Moscovich, Boaz Nadler, Clifford Spiegelman (2016). On the exact Berk-Jones statistics and their\n";
    cout << "        p-value calculation. Electronic Journal of Statistics. https://doi.org/10.1214/16-EJS1172\n";
    cout << "    [3] Amit Moscovich, Boaz Nadler (2017). Fast calculation of boundary crossing probabilities for Poisson processes.\n";
    cout << "        Statistics & Probability Letters. https://doi.org/10.1016/j.spl.2016.11.027\n";
}

double calculate_ecdf1_mns2016(int n, const vector<double>& b, const vector<double>& B)
{
    if ((b.size() > 0) && (B.size() == 0)) {
        return 1.0 - ecdf1_mns2016_upper(n, b);
    } else if ((b.size() == 0) && (B.size() > 0)) {
        return 1.0 - ecdf1_mns2016_lower(n, B);
    } else {
        print_usage();
        throw runtime_error("Expecting EITHER a lower or an upper boundary function when using the 'ecdf1-mns2016' command for computing a one-sided boundary crossing.\n");
    }
}

double calculate_ecdf1_new(int n, const vector<double>& b, const vector<double>& B)
{
    if ((b.size() > 0) && (B.size() == 0)) {
        return 1.0 - ecdf1_new_upper(n, b);
    } else if ((b.size() == 0) && (B.size() > 0)) {
        return 1.0 - ecdf1_new_lower(n, B);
    } else {
        print_usage();
        throw runtime_error("Expecting EITHER a lower or an upper boundary function when using the 'ecdf1-m2020' command for computing a one-sided boundary crossing.\n");
    }
}

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << '[';
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  return out;
}

double calculate_ecdf2_ks2001(int n, const vector<double>& b, const vector<double>& B)
{
    if ((b.size() == n) && (B.size() == n)) {
        return 1.0 - ecdf2(n, b, B, false);
    }

    if ((b.size() == 0) && (B.size() == n)) {
        std::vector<double> zeros_vector(n, 0.0);
        return 1.0 - ecdf2(n, zeros_vector, B, false);
    }

    if ((b.size() == n) && (B.size() == 0)) {
        std::vector<double> ones_vector(n, 1.0);
        return 1.0 - ecdf2(n, b, ones_vector, false);
    }

    throw runtime_error("Expecting either two boundary lists of length n or one list of length n and one of length zero");
}

double calculate_ecdf2_mn2017(int n, const vector<double>& b, const vector<double>& B)
{
    if ((b.size() == n) && (B.size() == n)) {
        return 1.0 - ecdf2(n, b, B, true);
    }

    if ((b.size() == 0) && (B.size() == n)) {
        std::vector<double> zeros_vector(n, 0.0);
        return 1.0 - ecdf2(n, zeros_vector, B, true);
    }

    if ((b.size() == n) && (B.size() == 0)) {
        std::vector<double> ones_vector(n, 1.0);
        return 1.0 - ecdf2(n, b, ones_vector, true);
    }
    throw runtime_error("Expecting either two boundary lists of length n or one list of length n and one of length zero");
}

static int handle_command_line_arguments(int argc, char* argv[])
{
    string command = string(argv[1]);

    long n = string_to_long(argv[2]);
    if (n < 0) {
        print_usage();
        throw runtime_error("n must be non-negative!");
    }

    string filename = string(argv[3]);
    pair<vector<double>, vector<double> > bounds = read_and_check_boundaries_file(filename);
    const vector<double>& b = bounds.first;
    const vector<double>& B = bounds.second;

    double result;
    if (command == "ecdf1-mns2016") {
        result = calculate_ecdf1_mns2016(n, b, B);
    } else if (command == "ecdf1-new") {
        result = calculate_ecdf1_new(n, b, B);
    } else if (command == "ecdf2-ks2001") {
        result = calculate_ecdf2_ks2001(n, b, B);
    } else if (command == "ecdf2-mn2017") {
        result = calculate_ecdf2_mn2017(n, b, B);
    } else {
        print_usage();
        throw runtime_error("Second command line argument must be one of: 'ecdf1-mns2016', 'ecdf1-new', 'ecdf2-ks2001', 'ecdf2-mn2017'.");
    }

    cout << result << endl;

    return 0;
}

int main(int argc, char* argv[])
{
    if (argc != 4) {
        print_usage();
        cout << "Error: Expecting 3 command line arguments!" << endl;
        return 1;
    }
    try {
        handle_command_line_arguments(argc, argv);
        return 0;
    } catch (ifstream::failure& e) {
        cout << "ifstream::failure exception caught:" << endl;
        cout << e.what() << endl;
        return 2;
    } catch (runtime_error& e) {
        cout << "Error:" << endl;
        cout << e.what() << endl;
        return 3;
    }
}
