#!/usr/bin/env python2.7
#
# This test script can be run directly from the shell, but using the "py.test" package gives nicer-looking output.

import subprocess

EPSILON = 0.01

def run(s):
    print 'Running %s' % s
    output = subprocess.check_output(s, shell=True)
    print 'Output: %s' % output.strip()
    return output

def test_ecdf():
    assert run('./bin/crossprob ecdf 2 tests/bounds2.txt').strip() == '0.25'
    assert run('./bin/crossprob ecdf 8 tests/bounds8.txt').strip() == '0.159471'
    assert run('./bin/crossprob ecdf 10 tests/bounds_cksplus_10.txt').strip() ==  '0.391076'
    assert run('./bin/crossprob ecdf 10 tests/bounds_cksminus_10.txt').strip() ==  '0.391076'
    
def test_poisson():
    assert run('./bin/crossprob poisson 2 tests/bounds2.txt').strip() == '0.661662'
    assert run('./bin/crossprob poisson 10 tests/bounds8.txt').strip() == '0.905357'
    assert run('./bin/crossprob poisson 7 tests/bounds_cksplus_10.txt').strip() ==  '0.231013'

def test_ecdf_one_sided():
    assert run('./bin/crossprob ecdf_one_sided 10 tests/bounds_cksminus_10.txt').strip() ==  '0.391076'
    assert run('./bin/crossprob ecdf_one_sided 10 tests/bounds_cksplus_10.txt').strip() ==  '0.391076'

def test_ecdf_one_sided_new():
    assert run('./bin/crossprob ecdf_one_sided_new 10 tests/bounds_cksminus_10.txt').strip() ==  '0.391076'
    assert run('./bin/crossprob ecdf_one_sided_new 10 tests/bounds_cksplus_10.txt').strip() ==  '0.391076'

def test_crossprob_mc_binomial():
    binomial_bounds2 = float(run('./bin/crossprob_mc ecdf 2 tests/bounds2.txt 1000000'))
    assert abs(binomial_bounds2 - 0.25) < EPSILON

    binomial_bounds8 = float(run('./bin/crossprob_mc ecdf 8 tests/bounds8.txt 1000000'))
    assert abs(binomial_bounds8 - 0.159471) < EPSILON

    binomial_cksplus_10 = float(run('./bin/crossprob_mc ecdf 10 tests/bounds_cksplus_10.txt 1000000'))
    assert abs(binomial_cksplus_10 - 0.391076) < EPSILON

    binomial_cksminus_10 = float(run('./bin/crossprob_mc ecdf 10 tests/bounds_cksminus_10.txt 1000000'))
    assert abs(binomial_cksminus_10 - 0.391076) < EPSILON

def test_crossprob_mc_poisson():
    poisson_bounds2 = float(run('./bin/crossprob_mc poisson 2 tests/bounds2.txt 1000000'))
    assert abs(poisson_bounds2 - 0.661662) < EPSILON

    poisson_bounds8 = float(run('./bin/crossprob_mc poisson 10 tests/bounds8.txt 1000000'))
    assert abs(poisson_bounds8 - 0.905357) < EPSILON

    poisson_cksplus_10 = float(run('./bin/crossprob poisson 7 tests/bounds_cksplus_10.txt'))
    assert abs(poisson_cksplus_10 - 0.231013) < EPSILON


def main():
    for (key,val) in globals().items():
        if key.startswith('test_'):
            val()

if __name__ == '__main__':
    main()

