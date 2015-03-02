#include <vector>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <cmath>
#include "one_sided_noncrossing_probability.hh"

using namespace std;

class Polynomial {
public:
    Polynomial(int max_degree);
    double get_coefficient(int degree) const;
    void set_coefficient(int degree, double value);
    double evaluate(double x) const;
    void integrate();
    friend ostream& operator<<(ostream& stream, const Polynomial& poly);
private:
    vector<double> coefficients;
    int degree;
};

Polynomial::Polynomial(int max_degree) :
    coefficients(max_degree+1, 0.0),
    degree(0)
{
    assert(max_degree >= 0);
}

double Polynomial::get_coefficient(int degree) const
{
    assert(degree >= 0);
    assert(degree <= this->degree);

    return coefficients[degree];

}

void Polynomial::set_coefficient(int degree, double value)
{
    assert(degree >= 0);
    assert(degree <= this->degree);

    coefficients[degree] = value;
}

double Polynomial::evaluate(double x) const
{
    double x_pow_i = 1.0;
    double result = 0.0;
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
        coefficients[i] = coefficients[i-1] / double(i);
    }
    coefficients[0] = 0.0;
    ++degree;
}

ostream& operator<<(ostream& stream, const Polynomial& poly)
{
    double coef;
    for (int i = poly.degree; i >= 1; --i) {
        coef = poly.coefficients[i];
        if (i != poly.degree) {
            stream << (coef >= 0 ? " + " : " - ");
        }
        stream << abs(coef) << " x";
        if (i > 1) {
            stream << "^" << i;
        }
    }
    coef = poly.coefficients[0];
    stream << (coef >= 0 ? " + " : " - ");
    stream << abs(coef);
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
        //cout << "i == " << i << ": " << p << endl;
        p.integrate();
        p.set_coefficient(0, -p.evaluate(upper_bound_steps[i]));
    }
    double integral_result = p.evaluate(1);
    return exp(lgamma(n+1) + log(integral_result));
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
