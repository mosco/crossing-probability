# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

"""
crossprob - crossing probabilities for one-dimensional empirical processes

Let X_1,...,X_n be a set of points sampled uniformly from the interval [0,1]
and let X_(1) <= X_(2) <= ... <= X_(n) be the sorted sample.

This module implements several algorithms for computing the probability
    Pr[for all i: b_i <= X_(i) <= B_i]        (1)

Equivalently, let the empirical cumulative distribution function of X_1,...,X_n be
    F_n(t) := \sum_i 1(X_i <= t).
The probability in Eq. (1) is equivalent to the probability that F_n(t) does not cross
a two-sided boundary,

The most versatile function is the two-sided crossing probability
    ecdf2(b, B, use_fft)
Its parameters are:
    n: the size of the sample
    b, B: two lists of length n of the boundaries in Eq. (1) above.
    use_fft: If true algorithm the O(n^2 logn) algorithm [MNS2016] is used,
             otherwise the O(n^3) algorithm of [KS2001]

Faster functions are available for the special case of a single boundary:
    ecdf1_new_b(b)
        Implements a new O(n^2) algorithm. B_i are implicitly assumed to be 1. 
    ecdf1_new_B(B)
        Implements the O(n^2) algorithm. b_i are implicitly assumed to be 0. 
    ecdf1_mns2016_b(b)
        Implements the O(n^2) algorithm of [MNS2016]. B_i are implicitly assumed to be 1. 
        Generally slower and less numerically stable than ecdf1_new_b()
    ecdf1_mns2016_B(B)
        Implements the O(n^2) algorithm of [MNS2016]. b_i are implicitly assumed to be 0. 
        Generally slower and less numerically stable than ecdf1_new_B()

EXAMPLES
    For a sample X_1, X_2, X_3 with order statistics X_(1) <= X_(2) <= X(3), the probability
        Pr[X_(1)<=0.7 and 0.15<=X_(2)<=0.9 and 0.5<=X_(3)<=0.8]
    may be computed using
        ecdf2([0, 0.15, 0.5], [0.7, 0.9, 0.8], True)

REFERENCES
    [KS2001] Estate Khmaladze, Eka Shinjikashvili (2001). Calculation of noncrossing probabilities for Poisson
             processes and its corollaries, Advances in Applied Probability. https://doi.org/10.1239/aap/1005091361
    [MNS2016] Amit Moscovich, Boaz Nadler, Clifford Spiegelman (2016). On the exact Berk-Jones statistics and their
              p-value calculation. Electronic Journal of Statistics. https://doi.org/10.1214/16-EJS1172
    [MN2017] Amit Moscovich, Boaz Nadler (2017). Fast calculation of boundary crossing probabilities for Poisson processes.
             Statistics & Probability Letters. https://doi.org/10.1016/j.spl.2016.11.027

MODULE REFERENCE
    https://github.com/mosco/crossing-probability

"""

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _crossprob
else:
    import _crossprob

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class SwigPyIterator(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _crossprob.delete_SwigPyIterator

    def value(self):
        return _crossprob.SwigPyIterator_value(self)

    def incr(self, n=1):
        return _crossprob.SwigPyIterator_incr(self, n)

    def decr(self, n=1):
        return _crossprob.SwigPyIterator_decr(self, n)

    def distance(self, x):
        return _crossprob.SwigPyIterator_distance(self, x)

    def equal(self, x):
        return _crossprob.SwigPyIterator_equal(self, x)

    def copy(self):
        return _crossprob.SwigPyIterator_copy(self)

    def next(self):
        return _crossprob.SwigPyIterator_next(self)

    def __next__(self):
        return _crossprob.SwigPyIterator___next__(self)

    def previous(self):
        return _crossprob.SwigPyIterator_previous(self)

    def advance(self, n):
        return _crossprob.SwigPyIterator_advance(self, n)

    def __eq__(self, x):
        return _crossprob.SwigPyIterator___eq__(self, x)

    def __ne__(self, x):
        return _crossprob.SwigPyIterator___ne__(self, x)

    def __iadd__(self, n):
        return _crossprob.SwigPyIterator___iadd__(self, n)

    def __isub__(self, n):
        return _crossprob.SwigPyIterator___isub__(self, n)

    def __add__(self, n):
        return _crossprob.SwigPyIterator___add__(self, n)

    def __sub__(self, *args):
        return _crossprob.SwigPyIterator___sub__(self, *args)
    def __iter__(self):
        return self

# Register SwigPyIterator in _crossprob:
_crossprob.SwigPyIterator_swigregister(SwigPyIterator)

class VectorDouble(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _crossprob.VectorDouble_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _crossprob.VectorDouble___nonzero__(self)

    def __bool__(self):
        return _crossprob.VectorDouble___bool__(self)

    def __len__(self):
        return _crossprob.VectorDouble___len__(self)

    def __getslice__(self, i, j):
        return _crossprob.VectorDouble___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _crossprob.VectorDouble___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _crossprob.VectorDouble___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _crossprob.VectorDouble___delitem__(self, *args)

    def __getitem__(self, *args):
        return _crossprob.VectorDouble___getitem__(self, *args)

    def __setitem__(self, *args):
        return _crossprob.VectorDouble___setitem__(self, *args)

    def pop(self):
        return _crossprob.VectorDouble_pop(self)

    def append(self, x):
        return _crossprob.VectorDouble_append(self, x)

    def empty(self):
        return _crossprob.VectorDouble_empty(self)

    def size(self):
        return _crossprob.VectorDouble_size(self)

    def swap(self, v):
        return _crossprob.VectorDouble_swap(self, v)

    def begin(self):
        return _crossprob.VectorDouble_begin(self)

    def end(self):
        return _crossprob.VectorDouble_end(self)

    def rbegin(self):
        return _crossprob.VectorDouble_rbegin(self)

    def rend(self):
        return _crossprob.VectorDouble_rend(self)

    def clear(self):
        return _crossprob.VectorDouble_clear(self)

    def get_allocator(self):
        return _crossprob.VectorDouble_get_allocator(self)

    def pop_back(self):
        return _crossprob.VectorDouble_pop_back(self)

    def erase(self, *args):
        return _crossprob.VectorDouble_erase(self, *args)

    def __init__(self, *args):
        _crossprob.VectorDouble_swiginit(self, _crossprob.new_VectorDouble(*args))

    def push_back(self, x):
        return _crossprob.VectorDouble_push_back(self, x)

    def front(self):
        return _crossprob.VectorDouble_front(self)

    def back(self):
        return _crossprob.VectorDouble_back(self)

    def assign(self, n, x):
        return _crossprob.VectorDouble_assign(self, n, x)

    def resize(self, *args):
        return _crossprob.VectorDouble_resize(self, *args)

    def insert(self, *args):
        return _crossprob.VectorDouble_insert(self, *args)

    def reserve(self, n):
        return _crossprob.VectorDouble_reserve(self, n)

    def capacity(self):
        return _crossprob.VectorDouble_capacity(self)
    __swig_destroy__ = _crossprob.delete_VectorDouble

# Register VectorDouble in _crossprob:
_crossprob.VectorDouble_swigregister(VectorDouble)


def ecdf2(b, B, use_fft):
    r"""ecdf2(VectorDouble b, VectorDouble B, bool use_fft) -> double"""
    return _crossprob.ecdf2(b, B, use_fft)

def ecdf1_mns2016_B(B):
    r"""ecdf1_mns2016_B(VectorDouble B) -> double"""
    return _crossprob.ecdf1_mns2016_B(B)

def ecdf1_mns2016_b(b):
    r"""ecdf1_mns2016_b(VectorDouble b) -> double"""
    return _crossprob.ecdf1_mns2016_b(b)

def ecdf1_new_B(B):
    r"""ecdf1_new_B(VectorDouble B) -> double"""
    return _crossprob.ecdf1_new_B(B)

def ecdf1_new_b(b):
    r"""ecdf1_new_b(VectorDouble b) -> double"""
    return _crossprob.ecdf1_new_b(b)


