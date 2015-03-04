#include <vector>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <sstream>

#include "one_sided_noncrossing_probability.hh"

using namespace std;

#define FLOAT_PRECISION_BITS 64

#if FLOAT_PRECISION_BITS == 128
    typedef __float128 FLOAT;
    #include <quadmath.h>
    #define EXP expq
    #define LOG logq
    #define GAMMA tgammaq
    #define LOG_GAMMA lgammaq
    #define POW powq
    #define FREXP frexpq
    #define LDEXP ldexpq
    #define MAX_EXP FLT128_MAX_EXP
    #define TO_STRING float128_to_string
    string float128_to_string(__float128 x)
    {
        static char s[1000];
        quadmath_snprintf(s, sizeof(s), "%.30Qg", x);
        return string(s);
    }
#elif FLOAT_PRECISION_BITS == 80
    typedef long double FLOAT;
    #include <cmath>
    #include <cfloat>
    #define EXP expl
    #define LOG logl
    #define GAMMA tgammal
    #define LOG_GAMMA lgammal
    #define POW powl
    #define FREXP frexpl
    #define LDEXP ldexpl
    #define MAX_EXP LDBL_MAX_EXP
    #define TO_STRING long_double_to_string
#elif FLOAT_PRECISION_BITS == 64
    typedef double FLOAT;
    #include <cmath>
    #include <cfloat>
    #define EXP exp
    #define LOG log
    #define GAMMA tgamma
    #define LOG_GAMMA lgamma
    #define POW pow
    #define FREXP frexp
    #define LDEXP ldexp
    #define MAX_EXP DBL_MAX_EXP
    #define TO_STRING double_to_string
#else
    #error FLOAT_PRECISION_BITS must be 64, 80 or 128.
#endif

string double_to_string(double x)
{
    stringstream ss;
    ss.precision(17); // http://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout
    ss << x;
    return ss.str();
}

string long_double_to_string(long double x)
{
    stringstream ss;
    ss.precision(22); // http://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout
    ss << x;
    return ss.str();
}

class Polynomial {
public:
    Polynomial(int max_degree);
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

void Polynomial::set_coefficient(int degree, FLOAT value)
{
    assert(degree >= 0);
    assert(degree <= this->degree);

    coefficients[degree] = value;
}

FLOAT Polynomial::evaluate(FLOAT x) const
{
    FLOAT result = 0.0;
    for (int i = 0; i < degree+1; ++i) {
        result += coefficients[i] * POW(x, i);
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

class PolynomialTranslatedMonomials {
public:
    PolynomialTranslatedMonomials(int max_degree);
    FLOAT get_multiplicative_coefficient(int degree) const;
    void set_multiplicative_coefficient(int degree, FLOAT multiplicative_coefficient);
    void set_additive_coefficient(int degree, FLOAT additive_coefficient);
    FLOAT evaluate(FLOAT x) const;
    void integrate();
    void ldexp_all_multiplicative_coefficients(int exp);
    friend ostream& operator<<(ostream& stream, const PolynomialTranslatedMonomials& poly);
    vector<FLOAT> multiplicative_coefficients;
    vector<FLOAT> additive_coefficients;
    int degree;
};

PolynomialTranslatedMonomials::PolynomialTranslatedMonomials(int max_degree) :
    multiplicative_coefficients(max_degree+1, 0.0),
    additive_coefficients(max_degree+1, 0.0),
    degree(0)
{
    assert(max_degree >= 0);
}

FLOAT PolynomialTranslatedMonomials::get_multiplicative_coefficient(int degree) const
{
    assert(degree >= 0);
    assert(degree <= this->degree);
    return multiplicative_coefficients[degree];
}

void PolynomialTranslatedMonomials::set_multiplicative_coefficient(int degree, FLOAT multiplicative_coefficient)
{
    assert(degree >= 0);
    assert(degree <= this->degree);

    multiplicative_coefficients[degree] = multiplicative_coefficient;
}

void PolynomialTranslatedMonomials::set_additive_coefficient(int degree, FLOAT additive_coefficient)
{
    assert(degree >= 0);
    assert(degree <= this->degree);

    additive_coefficients[degree] = additive_coefficient;
}

FLOAT PolynomialTranslatedMonomials::evaluate(FLOAT x) const
{
    FLOAT result = 0.0;
    for (int i = 0; i < degree+1; ++i) {
        result += multiplicative_coefficients[i] * POW(x+additive_coefficients[i], i);
    }
    return result;
}

void PolynomialTranslatedMonomials::integrate()
{
    assert(degree <= (int)multiplicative_coefficients.size());
    for (int i = degree+1; i >= 1; --i) {
        multiplicative_coefficients[i] = multiplicative_coefficients[i-1] / ((FLOAT)i);
        additive_coefficients[i] = additive_coefficients[i-1];
    }
    additive_coefficients[1] = 0.0;
    additive_coefficients[0] = 0.0;
    multiplicative_coefficients[0] = 0.0;
    ++degree;
}

void PolynomialTranslatedMonomials::ldexp_all_multiplicative_coefficients(int exp)
{
    for (int i = 0; i < degree+1; ++i) {
        multiplicative_coefficients[i] = LDEXP(multiplicative_coefficients[i], exp);
    }
}

static void print_FLOAT_array(const FLOAT* arr, int n)
{
    for (int i = 0; i < n; ++i) {
        cout << TO_STRING(arr[i]) << ", ";
    }
    cout << endl;
}

ostream& operator<<(ostream& stream, const Polynomial& poly)
{
    for (int i = poly.degree; i >= 0; --i) {
        FLOAT coef = poly.coefficients[i];
        if (coef >= 0) {
            stream << (i != poly.degree ? " + " : "") << TO_STRING(coef);
            
        } else {
            stream << (i != poly.degree ? " - " : "") << TO_STRING(-coef);
        }
        if (i >= 1) {
            stream << " x^" << i;
        }
    }
    return stream;
}
ostream& operator<<(ostream& stream, const PolynomialTranslatedMonomials& poly)
{
    for (int i = poly.degree; i >= 0; --i) {
        FLOAT coef = poly.multiplicative_coefficients[i];
        if (coef >= 0) {
            stream << (i != poly.degree ? " + " : "") << TO_STRING(coef);
            
        } else {
            stream << (i != poly.degree ? " - " : "") << TO_STRING(-coef);
        }
        if (i >= 1) {
            stream << " (x + " << TO_STRING(poly.additive_coefficients[i]) << ")^" << i;
        }
    }
    return stream;
}

double binomial_process_upper_noncrossing_probability_naive(int n, const vector<double>& upper_bound_steps)
{
    if ((int)upper_bound_steps.size() < n) {
        throw runtime_error("Binomial process b(t) must cross upper boundary h(t) since h(1) < n and b(t) = n");
    }
    Polynomial p(n);
    p.set_coefficient(0, 1.0);
    for (int i = 0; i < n; ++i) {
        cout << "============= i == " << i << ": " << p << endl;
        //cout << "Before integration: ";
        //print_FLOAT_array(&p.coefficients[0], p.degree+1);
        p.integrate();
        //cout << "After integration: ";
        //print_FLOAT_array(&p.coefficients[0], p.degree+1);
        //cout << "After integration (<<): " << p << endl;
        //cout << "Setting constant to -p.evaluate(" << (double)upper_bound_steps[i] << ") == " << (double)-p.evaluate(upper_bound_steps[i]) << endl;
        p.set_coefficient(0, -p.evaluate(upper_bound_steps[i]));
    }
    FLOAT integral_result = p.evaluate(1);
    cout << "Integral result: " << TO_STRING(integral_result) << endl;
    cout << "LOG_GAMMA(n+1) == " << TO_STRING(LOG_GAMMA(n+1)) << endl;
    cout << "LOG(integral_result) == " << TO_STRING(LOG(integral_result)) << endl;
    cout << "EXP(LOG_GAMMA(n+1) + LOG(integral_result))" << TO_STRING(EXP(LOG_GAMMA(n+1) + LOG(integral_result))) << endl;
    return EXP(LOG_GAMMA(n+1) + LOG(integral_result));
}

double binomial_process_upper_noncrossing_probability_lazy_integration(int n, const vector<double>& upper_bound_steps)
{
    if ((int)upper_bound_steps.size() < n) {
        throw runtime_error("Binomial process b(t) must cross upper boundary h(t) since h(1) < n and b(t) = n");
    }
    vector<FLOAT> constant_terms(n+1, 0);
    constant_terms[0] = 1.0;
    for (int i = 0; i < n; ++i) {
        FLOAT new_constant_term = 0;
        for (int j = 0; j <=i; j++) {
            new_constant_term -= constant_terms[j] * POW(upper_bound_steps[i], i-j+1) / GAMMA(i-j+2); // Li^{i-j+1}/(i-j+1)!
        }
        constant_terms[i+1] = new_constant_term;
    }
    cout << "constant_terms:\n";
    print_FLOAT_array(&constant_terms[0], n+1);

    FLOAT evaluation_at_1 = 0.0;
    for (int i = 0; i < n+1; ++i) {
        evaluation_at_1 += constant_terms[i] / GAMMA(n+1-i);
    }
    FLOAT integral_result = evaluation_at_1;
    cout << "Integral result: " << TO_STRING(integral_result) << endl;
    cout << "LOG_GAMMA(n+1) == " << TO_STRING(LOG_GAMMA(n+1)) << endl;
    cout << "LOG(integral_result) == " << TO_STRING(LOG(integral_result)) << endl;
    cout << "EXP(LOG_GAMMA(n+1) + LOG(integral_result))" << TO_STRING(EXP(LOG_GAMMA(n+1) + LOG(integral_result))) << endl;
    return EXP(LOG_GAMMA(n+1) + LOG(integral_result));
}

double binomial_process_upper_noncrossing_probability_trans(int n, const vector<double>& upper_bound_steps)
{
    cout << "Using translated polynomials! precision: " << FLOAT_PRECISION_BITS << endl;
    if ((int)upper_bound_steps.size() < n) {
        throw runtime_error("Binomial process b(t) must cross upper boundary h(t) since h(1) < n and b(t) = n");
    }
    PolynomialTranslatedMonomials p(n+1);
    p.set_multiplicative_coefficient(0, 1.0);

    for (int i = 0; i < (int)upper_bound_steps.size(); ++i) {
        cout << "============= i == " << i << ": " << p << endl;
        p.integrate();
        p.set_additive_coefficient(1, -upper_bound_steps[i]);
        p.set_multiplicative_coefficient(0, -p.evaluate(upper_bound_steps[i]));
    }

    FLOAT integral_result = p.evaluate(1);
    cout << "Integral result: " << TO_STRING(integral_result) << endl;
    cout << "LOG_GAMMA(n+1) == " << TO_STRING(LOG_GAMMA(n+1)) << endl;
    cout << "LOG(integral_result) == " << TO_STRING(LOG(integral_result)) << endl;
    cout << "EXP(LOG_GAMMA(n+1) + LOG(integral_result))" << TO_STRING(EXP(LOG_GAMMA(n+1) + LOG(integral_result))) << endl;
    return EXP(LOG_GAMMA(n+1) + LOG(integral_result));
}

FLOAT maximum_multiplicative_coefficient(const PolynomialTranslatedMonomials& p)
{
    assert(p.degree >= 1);
    FLOAT m = p.get_multiplicative_coefficient(0);
    for (int i = 0; i < p.degree+1; ++i) {
        FLOAT x = p.get_multiplicative_coefficient(i);
        if (x > m) {
            m = x;
        }
    }
    return m;
}

int extract_exponent(FLOAT x)
{
    static int exponent;
    FREXP(x, &exponent);
    return exponent;
}

double binomial_process_upper_noncrossing_probability(int n, const vector<double>& upper_bound_steps)
{
    cout << "Using translated polynomials! precision: " << FLOAT_PRECISION_BITS << endl;
    if ((int)upper_bound_steps.size() < n) {
        throw runtime_error("Binomial process b(t) must cross upper boundary h(t) since h(1) < n and b(t) = n");
    }
    PolynomialTranslatedMonomials p(n+1);
    p.set_multiplicative_coefficient(0, 1.0);

    int total_exponent_delta = 0;
    for (int i = 0; i < (int)upper_bound_steps.size(); ++i) {
        //cout << "============= i == " << i << ": " << p << endl;
        p.integrate();
        p.set_additive_coefficient(1, -upper_bound_steps[i]);
        p.set_multiplicative_coefficient(0, -p.evaluate(upper_bound_steps[i]));

        //cout << "Before expfix: " << p << endl;

        FLOAT max_coef = maximum_multiplicative_coefficient(p);
        int highest_exponent = extract_exponent(max_coef);
        if (highest_exponent < 0) {
            //cout << "highest exponent: " << highest_exponent << endl;
            p.ldexp_all_multiplicative_coefficients(0 - highest_exponent);
            total_exponent_delta += 0 - highest_exponent;
        }
    }
    FLOAT integral_result = p.evaluate(1);
    FLOAT log_integral_result = LOG(integral_result);
    cout << "Integral result: (no exponent fix) " << TO_STRING(integral_result) << endl;
    cout << "Exponent delta: " << total_exponent_delta << endl;
    cout << "LOG_GAMMA(n+1) == " << TO_STRING(LOG_GAMMA(n+1)) << endl;
    cout << "log_integral_result (no exponent fix) == " << TO_STRING(log_integral_result) << endl;
    FLOAT noncrossing_probability = EXP(LOG_GAMMA(n+1) + log_integral_result - total_exponent_delta*LOG(2));
    cout << "Final result: " << TO_STRING(noncrossing_probability) << endl;
    return noncrossing_probability;
}

//Plain polynomials:
//1
//x-L0
//1/2 x^2 - L0 x - (1/2 L1^2 - L0L1)
//1/6 x^3 - 1/2 L0 x^2 - (1/2 L1^2 - L0L1)x - 1/6 L2^3 + 1/2 L0 L2^2 + 1/2 L1^2L2 - L0L1L2
//
//
//Translated polynomials:
//1
//(x-L0) 
//1/2 (x-L0)^2 + 0*(x-L1) - (1/2 (L1-L0)^2) = 1/2 x^2 - L0 x - 1/2 L1^2  + L0L1



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
