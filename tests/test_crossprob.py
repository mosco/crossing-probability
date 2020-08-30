#!/usr/bin/env python3
#
# This test script can be run directly from the shell, but using the "py.test" package gives nicer-looking output.

import subprocess

EPSILON = 0.01

def run(s):
    print('Running ', s)
    output = subprocess.check_output(s, shell=True)
    print('Output:', output.strip())
    return output

def test_ecdf2ks2001():
    assert run('./bin/crossprob ecdf2-ks2001 tests/bounds_0_1.txt').strip() == b'1'
    assert run('./bin/crossprob ecdf2-ks2001 tests/bounds2.txt').strip() == b'0.75'
    assert run('./bin/crossprob ecdf2-ks2001 tests/bounds8.txt').strip() == b'0.840529'
    assert run('./bin/crossprob ecdf2-ks2001 tests/bounds_cksplus_10.txt').strip() ==  b'0.608924'
    assert run('./bin/crossprob ecdf2-ks2001 tests/bounds_cksminus_10.txt').strip() ==  b'0.608924'
    
#def test_poisson():
#    assert run('./bin/crossprob poisson 2 tests/bounds2.txt').strip() == b'0.661662'
#    assert run('./bin/crossprob poisson 10 tests/bounds8.txt').strip() == b'0.905357'
#    assert run('./bin/crossprob poisson 7 tests/bounds_cksplus_10.txt').strip() ==  b'0.231013'

def test_ecdf1mns2016():
    assert run('./bin/crossprob ecdf1-mns2016 tests/bounds_0.txt').strip() ==  b'1'
    assert run('./bin/crossprob ecdf1-mns2016 tests/bounds__1.txt').strip() ==  b'1'
    assert run('./bin/crossprob ecdf1-mns2016 tests/bounds_cksminus_10.txt').strip() ==  b'0.608924'
    assert run('./bin/crossprob ecdf1-mns2016 tests/bounds_cksplus_10.txt').strip() ==  b'0.608924'

def test_ecdf2mn2017():
    assert run('./bin/crossprob ecdf2-mn2017 tests/bounds_0_1.txt').strip() == b'1'
    assert run('./bin/crossprob ecdf2-mn2017 tests/bounds2.txt').strip() == b'0.75'
    assert run('./bin/crossprob ecdf2-mn2017 tests/bounds8.txt').strip() == b'0.840529'
    assert run('./bin/crossprob ecdf2-mn2017 tests/bounds_cksminus_10.txt').strip() ==  b'0.608924'
    assert run('./bin/crossprob ecdf2-mn2017 tests/bounds_cksplus_10.txt').strip() ==  b'0.608924'

def test_ecdf1m2020():
    assert run('./bin/crossprob ecdf1-new tests/bounds_0.txt').strip() ==  b'1'
    assert run('./bin/crossprob ecdf1-new tests/bounds__1.txt').strip() ==  b'1'
    assert run('./bin/crossprob ecdf1-new tests/bounds_cksminus_10.txt').strip() ==  b'0.608924'
    assert run('./bin/crossprob ecdf1-new tests/bounds_cksplus_10.txt').strip() ==  b'0.608924'

def test_crossprob_mc_binomial():
    binomial_bounds_0 = float(run('./bin/crossprob_mc ecdf tests/bounds_0.txt 1000'))
    assert binomial_bounds_0 == 1
    binomial_bounds__1 = float(run('./bin/crossprob_mc ecdf tests/bounds__1.txt 1000'))
    assert binomial_bounds__1 == 1

    binomial_bounds2 = float(run('./bin/crossprob_mc ecdf tests/bounds2.txt 1000000'))
    assert abs(binomial_bounds2 - 0.75) < EPSILON

    binomial_bounds8 = float(run('./bin/crossprob_mc ecdf tests/bounds8.txt 1000000'))
    assert abs(binomial_bounds8 - 0.840529) < EPSILON

    binomial_cksplus_10 = float(run('./bin/crossprob_mc ecdf tests/bounds_cksplus_10.txt 1000000'))
    assert abs(binomial_cksplus_10 - 0.608924) < EPSILON

    binomial_cksminus_10 = float(run('./bin/crossprob_mc ecdf tests/bounds_cksminus_10.txt 1000000'))
    assert abs(binomial_cksminus_10 - 0.608924) < EPSILON

#def test_crossprob_mc_poisson():
#    poisson_bounds2 = float(run('./bin/crossprob_mc poisson 2 tests/bounds2.txt 1000000'))
#    assert abs(poisson_bounds2 - 0.661662) < EPSILON
#
#    poisson_bounds8 = float(run('./bin/crossprob_mc poisson 10 tests/bounds8.txt 1000000'))
#    assert abs(poisson_bounds8 - 0.905357) < EPSILON
#
#    poisson_cksplus_10 = float(run('./bin/crossprob poisson 7 tests/bounds_cksplus_10.txt'))
#    assert abs(poisson_cksplus_10 - 0.231013) < EPSILON


def main():
    for (key,val) in globals().items():
        if key.startswith('test_'):
            val()

if __name__ == '__main__':
    main()

