# Benchmarks for testing different methods for one-sided boundary crossing of the empirical CDF.
#
# The plot_all function generates the figures from the paper:
#     "Fast calculation of p-values for one-sided Kolmogorov-Smirnov type statistics"
#     https://arxiv.org/abs/2009.04954
#
# To run these benchmarks you will need to have SciPy and Matplotlib installed, then follow the instructions in the README.txt (located in the main dir) to build the C++ code and Python module.
#
# Amit Moscovich, Tel Aviv University, 2022.

import sys
import time
from numpy import arange, array, asarray, double, longdouble, absolute, inf
from scipy.special import betainc, betaincinv

import matplotlib
import matplotlib.pyplot as plt

import cpuinfo # pip install py-cpuinfo


import crossprob
from utils import Timer
import pickler
from latex_fig_utils import RCPARAMS_LATEX_SINGLE_COLUMN_LARGE, RCPARAMS_LATEX_SINGLE_COLUMN_LARGE_TALL, RCPARAMS_LATEX_SINGLE_COLUMN_WIDE, save_figure

def _uniform_order_statistic_cdf(i, n, x):
    """
    _uniform_order_statistic_cdf(i, n, x) -> Pr[U_(i) < x]
    
    Let U_1, ..., U_n ~ Uniform[0,1] be n independent random variables
    and let U_(1) < ... < U_(n) denote the same variables in sorted order.  Then
        U_{(i)} ~ Beta(i, n-i+1)

    This function returns the Cumulative Distribution function of U_(i),
    i.e. the return value is Pr[U_(i) < x]

    This function also works for numpy array inputs"""
    return betainc(i, n-i+1, x)


def Mn_plus_statistic(uniform01_samples):
    assert all(0 <= x <= 1 for x in uniform01_samples)
    X = asarray(uniform01_samples)
    n = len(uniform01_samples)
    X.sort()
    p_values = _uniform_order_statistic_cdf(arange(1,n+1), n, X)
    return p_values.min()


def Mn_plus_bounds(n, Mn_plus_value):
    indices = arange(1, n+1)
    return betaincinv(indices, n-indices+1, Mn_plus_value)


def Mn_plus_distribution(n, x):
    """Probability of having Mn_plus < x under the null hypothesis that X_1,...X_n ~ U[0,1]"""
    b_bounds = Mn_plus_bounds(n, x)
    return 1.0 - crossprob.ecdf1_new_b(b_bounds)
    #return 1.0 - crossprob.ecdf2(b_bounds, [1]*len(b_bounds), True)

N_BINARY_SEARCH_STEPS_MAX = 110 
EPSILON = 1e-11

def inverse_Mn_plus(n, alpha, debug_prints):
    """
    Binary search for x such that
        Mn_plus_distribution(n,x) = alpha
    """
    low = 0.0
    high = 1.0 
    last_range = None
    for i in range(N_BINARY_SEARCH_STEPS_MAX):
        if last_range == (low, high):
            if debug_prints:
                print('last_range == (low, high)')
            break
        last_range = (low, high)

        mid = (low+high)/2.0
        if debug_prints:
            print(f'{i}: low={low:.25f} mid={mid:.25f} high={high:.25f}')

        mid_alpha = Mn_plus_distribution(n, mid)
        if debug_prints:
            print(f'Pr[M_{n} < {mid:.25f}] = {mid_alpha:.25f}')
        if mid_alpha < alpha:
            low = mid
        else:
            high = mid

        if low >= high:
            if debug_prints:
                print('low >= high')
            break
    assert i < N_BINARY_SEARCH_STEPS_MAX-1, "Reached N_BINARY_SEARCH_STEPS_MAX => did not converge!"
    result = mid
    alpha_bounds = Mn_plus_bounds(n, result)
    nocross_probability = 1.0 - crossprob.ecdf1_new_b(alpha_bounds)
    relative_error = absolute(alpha - nocross_probability)/alpha 
    if relative_error > EPSILON:
        print('Large relative error!', relative_error)

    return result

class Timer(object):
    def __init__(self, text=None):
        if text is not None:
            print('%s: ' % text, end='') 
            sys.stdout.flush()
        self.text = text
        self.start_clock = time.clock()
        self.start_time = time.time()

    def stop(self):
        self.stop_time = time.time()
        self.elapsed = self.stop_time - self.start_time

        self.stop_clock = time.clock()
        self.elapsed_cpu_time = self.stop_clock - self.start_clock

    def print_elapsed(self):
        print('Wall time: %.3f seconds.  CPU time: %.3f seconds.' % (self.elapsed, self.elapsed_cpu_time))

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        self.stop()
        if self.text is not None:
            self.print_elapsed()


def precalc_Mn_plus_alpha_level_thresholds(alpha, n_range):
    thresholds = []
    N_RANGE = list(n_range)
    for n in N_RANGE:
        print('n:', n)
        with Timer():
            threshold = inverse_Mn_plus(n, alpha, debug_prints=False)
        thresholds.append((n, threshold))
        pickler.dump(f'Mn_plus_thresholds_{alpha}', n_range=N_RANGE, thresholds=thresholds)


def precalc_timings_Mn_plus(alpha, max_ks2001_n, max_mns2016_n, num_reps):
    N_TIMING_REPETITIONS = 3
    CPU_BRAND_STRING = cpuinfo.get_cpu_info()['brand_raw']
    print(f'Running on {CPU_BRAND_STRING}')
    
    d = pickler.load(f'Mn_plus_thresholds_{alpha}')
    
    timings_ks2001 = [inf]*len(d.n_range)
    results_ks2001 = [None]*len(d.n_range)

    timings_mns2016 = [inf]*len(d.n_range)
    results_mns2016 = [None]*len(d.n_range)


    timings_mn2017 = [inf]*len(d.n_range)
    results_mn2017 = [None]*len(d.n_range)

    timings_new = [inf]*len(d.n_range)
    results_new = [None]*len(d.n_range)

    for rep in range(num_reps):
        print(f'==== Run {rep} =======================================')
        for (i, (n, threshold)) in enumerate(d.thresholds):
            print(f'---- n={n} -----------------')
            with Timer('Computing bounds'):
                bounds = Mn_plus_bounds(n, threshold)

            if n <= max_ks2001_n:
                with Timer('KS2001') as t:
                    results_ks2001[i] = 1.0-crossprob.ecdf2(bounds, [1.0]*n, False)
                timings_ks2001[i] = min(timings_ks2001[i], t.elapsed)

            if n <= max_mns2016_n:
                with Timer('MNS2016') as t:
                    results_mns2016[i] = 1.0-crossprob.ecdf1_mns2016_b(bounds)
                timings_mns2016[i] = min(timings_mns2016[i], t.elapsed)

            with Timer('MN2017') as t:
                results_mn2017[i] = 1.0-crossprob.ecdf2(bounds, [1.0]*n, True)
            timings_mn2017[i] = min(timings_mn2017[i], t.elapsed)

            with Timer('NEW') as t:
                results_new[i] = 1.0-crossprob.ecdf1_new_b(bounds)
            timings_new[i] = min(timings_new[i], t.elapsed)

            pickler.dump(f'Mn_plus_timings_{alpha}',
                n_range=d.n_range,
                cpu = CPU_BRAND_STRING,
                num_reps = num_reps,
                ks2001 = (timings_ks2001, results_ks2001),
                mns2016 = (timings_mns2016, results_mns2016),
                mn2017 = (timings_mn2017, results_mn2017),
                new = (timings_new, results_new))


def precalc_timings_Mn_plus_largescale(threshold, n_range):
    CPU_BRAND_STRING = cpuinfo.get_cpu_info()['brand_raw']
    print(f'Running on {CPU_BRAND_STRING}')
    
    timings_mn2017 = [inf]*len(n_range)
    results_mn2017 = [None]*len(n_range)

    timings_new = [inf]*len(n_range)
    results_new = [None]*len(n_range)

    for (i,n) in enumerate(n_range):
        print(f'---- n={n} -----------------')

        with Timer('Computing bounds'):
            bounds = Mn_plus_bounds(n, threshold)

        with Timer('NEW') as t:
            results_new[i] = 1.0-crossprob.ecdf1_new_b(bounds)
        timings_new[i] = min(timings_new[i], t.elapsed)

        with Timer('MN2017') as t:
            results_mn2017[i] = 1.0-crossprob.ecdf2(bounds, [1.0]*n, True)
        timings_mn2017[i] = min(timings_mn2017[i], t.elapsed)


        pickler.dump(f'Mn_plus_timings_largescale',
            n_range=n_range,
            threshold=threshold,
            cpu = CPU_BRAND_STRING,
            mn2017 = (timings_mn2017, results_mn2017),
            new = (timings_new, results_new))


def plot_timings_Mn_plus(alpha, xlim=None, ylim=None, xticks=None, yticks=None, name=None):
    d = pickler.load(f'Mn_plus_timings_{alpha}')

    with matplotlib.rc_context(rc = RCPARAMS_LATEX_SINGLE_COLUMN_LARGE_TALL):
        fig, ax = plt.subplots()

        ax.set_yscale('log')
        ax.plot(d.n_range, d.ks2001[0], label='KS (2001)')

        mns2016 = [(n,t) for (n,t) in zip(d.n_range, d.mns2016[0]) if n <= 50000]
        ax.plot(list(zip(*mns2016))[0], list(zip(*mns2016))[1], label='MNS (2016)')

        ax.plot(d.n_range, d.mn2017[0], label='MN (2017)')

        ax.plot(d.n_range, d.new[0], label='New')
        ax.legend(loc='lower right')
        ax.grid(which='both')
        #ax.set_title('Benchmarks')

        ax.minorticks_on()

        ax.set_ylabel('Runtime [sec]')
        ax.set_xlabel('n')


        ax.set_xticks(xticks)
        ax.set_xticks(xticks, minor=True)
        ax.set_xticklabels([f'{n:,}' for n in xticks], rotation=90)

        if yticks is not None:
            ax.set_yticks(yticks)

        # Slanted ticks at the bottom and top
        ax.tick_params(axis="x", bottom=True, top=True, labelbottom=True, labeltop=True)
        plt.setp([tick.label1 for tick in ax.xaxis.get_major_ticks()], rotation=45, ha="right", va="center", rotation_mode="anchor")
        plt.setp([tick.label2 for tick in ax.xaxis.get_major_ticks()], rotation=45, ha="left", va="center",rotation_mode="anchor")

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if name is None:
            name = f'Mn_plus_timings_{alpha}'

        save_figure(fig, name)


def plot_timings_Mn_plus_largescale(xlim=None, ylim=None, xticks=None, yticks=None, name=None):
    d = pickler.load(f'Mn_plus_timings_largescale')

    with matplotlib.rc_context(rc = RCPARAMS_LATEX_SINGLE_COLUMN_LARGE):
        fig, ax = plt.subplots()

        ax.set_yscale('log')

        # Advance axes.prop_cycle twice, so that this plot is color and marker-compatible with plots generated by plot_timings_Mn_plus
        ax.plot([], [])
        ax.plot([], [])
        ax.plot(d.n_range, d.mn2017[0], label='MN (2017)')
        ax.plot(d.n_range, d.new[0], label='New')
        ax.legend(loc='lower right')
        ax.grid(which='both')
        #ax.set_title('Benchmarks')

        ax.minorticks_on()

        ax.set_ylabel('Runtime [sec]')
        #ax.set_xlabel('n')


        if xticks is not None:
            ax.set_xticks(xticks)
            ax.set_xticks(xticks, minor=True)
            ax.set_xticklabels([f'{n:,}' for n in xticks], rotation=90)

        if yticks is not None:
            ax.set_yticks(yticks)

        # Slanted ticks at the bottom and top
        ax.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True)
        plt.setp([tick.label1 for tick in ax.xaxis.get_major_ticks()], rotation=45, ha="right", va="center", rotation_mode="anchor")
        plt.setp([tick.label2 for tick in ax.xaxis.get_major_ticks()], rotation=45, ha="left", va="center",rotation_mode="anchor")

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if name is None:
            name = f'Mn_plus_timings_largescale'

        save_figure(fig, name)

def plot_relative_errors():
    d = pickler.load(f'Mn_plus_timings_largescale')

    with matplotlib.rc_context(rc = RCPARAMS_LATEX_SINGLE_COLUMN_WIDE):
        fig, ax = plt.subplots()

        #ax.set_yscale('log')

        res_mn2017 = array(d.mn2017[1])
        res_new = array(d.new[1])
        relative_error = absolute(res_new - res_mn2017) / res_mn2017

        ax.plot(d.n_range, relative_error, '.-k', label='Relative error')
        #ax.legend(loc='lower right')
        ax.grid(which='both')
        #ax.set_title('Benchmarks')

        #ax.minorticks_on()

        ax.set_ylabel('Relative error')
        ax.set_xlabel('n')

        xticks = d.n_range
        ax.set_xticks(xticks)
        ax.set_xticks(xticks, minor=True)
        ax.set_xticklabels([f'{n:,}' for n in xticks], rotation=90)

        ax.set_xlim([d.n_range[0], d.n_range[-1]])

        ax.set_ylim([0, 1.2e-9])
        #ax.set_yticks(arange(0, 1.5e-9, 0.25e-9))
        ax.set_yticks(arange(0, 1.4e-9, 0.2e-9))

        # Slanted ticks at the bottom and top
        ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
        plt.setp([tick.label1 for tick in ax.xaxis.get_major_ticks()], rotation=45, ha="right", va="center", rotation_mode="anchor")
        #plt.setp([tick.label2 for tick in ax.xaxis.get_major_ticks()], rotation=45, ha="left", va="center",rotation_mode="anchor")

        save_figure(fig, 'Mn_plus_largescale_relative_error')


ALPHA = 0.05
N_RANGE = list(range(5000,105000,5000))

# KS2001 is too slow for testing large sample sizes.
KS2001_MAX_N = 35000

# MNS2016 starts producing inaccurate results for n>30,000 and also produces lots of NaN values, which slow down the computation.
MNS2016_MAX_N = 30000
NUM_REPS = 3

LARGESCALE_THRESHOLD = 0.0005437813829870565 # 0.05-level threshold for M_n_plus with n=100,000
LARGESCALE_RANGE = range(100000,1100000,100000)

def precalc_all():
    precalc_Mn_plus_alpha_level_thresholds(ALPHA, N_RANGE)
    precalc_timings_Mn_plus(ALPHA, KS2001_MAX_N, MNS2016_MAX_N, NUM_REPS) 

    precalc_timings_Mn_plus_largescale(LARGESCALE_THRESHOLD, LARGESCALE_RANGE)


def plot_all():
    plot_timings_Mn_plus(ALPHA, xlim=(5000, N_RANGE[-1]), ylim=(0.01,1000), xticks=range(5000,105000,5000), yticks=[0.01, 0.1, 1, 10, 100, 1000], name=f'Mn_plus_timings_{ALPHA}')
    plot_timings_Mn_plus_largescale(xlim=(100000,1000000), ylim=(10, 100000), xticks=LARGESCALE_RANGE)
    plot_relative_errors()
    

