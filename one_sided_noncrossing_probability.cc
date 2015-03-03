#include <vector>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <quadmath.h>
#include <math.h>
#include "one_sided_noncrossing_probability.hh"

using namespace std;

typedef __float128 FLOAT;
#define EXP expq
#define LOG logq
#define LOG_GAMMA lgammaq

//typedef double FLOAT;
//#define EXP exp
//#define LOG log
//#define LOG_GAMMA lgamma

class Polynomial {
public:
    Polynomial(int max_degree);
    FLOAT get_coefficient(int degree) const;
    void set_coefficient(int degree, FLOAT value);
    FLOAT evaluate(FLOAT x) const;
    void integrate();
    friend ostream& operator<<(ostream& stream, const Polynomial& poly);
    vector<FLOAT> coefficients;
    int degree;
};

Polynomial::Polynomial(int max_degree) :
    coefficients(max_degree+1, 0.0),
    degree(0)
{
    assert(max_degree >= 0);
}

FLOAT Polynomial::get_coefficient(int degree) const
{
    assert(degree >= 0);
    assert(degree <= this->degree);

    return coefficients[degree];

}

void Polynomial::set_coefficient(int degree, FLOAT value)
{
    assert(degree >= 0);
    assert(degree <= this->degree);

    coefficients[degree] = value;
}

FLOAT Polynomial::evaluate(FLOAT x) const
{
    FLOAT x_pow_i = 1.0;
    FLOAT result = 0.0;
    for (int i = 0; i < degree+1; ++i) {
        result += coefficients[i] * x_pow_i;
        x_pow_i *= x;
    }
    return result;
}

void Polynomial::integrate()
{
    assert(degree <= (int)coefficients.size());
    for (int i = degree+1; i >= 1; --i) {
        coefficients[i] = coefficients[i-1] / ((FLOAT)i);
    }
    coefficients[0] = 0.0;
    ++degree;
}

//static void print_FLOAT_array(const FLOAT* arr, int n)
//{
//    for (int i = 0; i < n; ++i) {
//        cout << arr[i] << ", ";
//    }
//    cout << endl;
//}

ostream& operator<<(ostream& stream, const Polynomial& poly)
{
    for (int i = poly.degree; i >= 0; --i) {
        FLOAT coef = poly.coefficients[i];
        if (coef >= 0) {
            stream << (i != poly.degree ? " + " : "") << (double)coef;
        } else {
            stream << (i != poly.degree ? " - " : "") << (double)-coef;
        }
        if (i >= 1) {
            stream << " x^" << i;
        }
    }
    return stream;
}

double binomial_process_upper_noncrossing_probability(int n, const vector<double>& upper_bound_steps)
{
    if ((int)upper_bound_steps.size() < n) {
        throw runtime_error("Binomial process b(t) must cross upper boundary h(t) since h(1) < n and b(t) = n");
        return 0;
    }
    Polynomial p(n);
    p.set_coefficient(0, 1.0);
    for (int i = 0; i < (int)upper_bound_steps.size(); ++i) {
        cout << "i == " << i << ": " << p << endl;
        //cout << "Before integration: ";
        //print_FLOAT_array(&p.coefficients[0], p.degree+1);
        p.integrate();
        //cout << "After integration: ";
        //print_FLOAT_array(&p.coefficients[0], p.degree+1);
        //cout << "After integration (<<): " << p << endl;
        p.set_coefficient(0, -p.evaluate(upper_bound_steps[i]));
    }
    FLOAT integral_result = p.evaluate(1);
    return EXP(LOG_GAMMA(n+1) + LOG(integral_result));
}

double binomial_process_lower_noncrossing_probability(int n, const vector<double>& lower_bound_steps)
{
    if ((int)lower_bound_steps.size() > n) {
        throw runtime_error("Binomial process b(t) must cross lower boundary g(t) since g(1) > n and b(t) = n.");
        return 0;
    }

    vector<double> symmetric_steps(n, 0.0);
    for (int i = n-lower_bound_steps.size(); i < n; ++i) {
        symmetric_steps[i] = 1.0 - lower_bound_steps[lower_bound_steps.size() - 1 - i];
    }

    return binomial_process_upper_noncrossing_probability(n, symmetric_steps);
}
