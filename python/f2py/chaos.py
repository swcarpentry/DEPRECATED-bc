#!python
# chaos.py
import pylab as pl
import numpy as np

# we import the fortran extension module here
import _chaos

# here is the logistic function
# this uses some advanced Python features.
# Logistic is a function that returns another function.
# This is known as a 'closure' and is a very powerful feature.
def logistic(r):
    def _inner(x):
        return r * x * (1.0 - x)
    return _inner

def sine(r):
    from math import sin, pi
    def _inner(x):
        return r * sin(pi * x)
    return _inner

def driver(func, lower, upper, N=400):
    # X will scan over the parameter value.
    X = np.linspace(lower, upper, N)
    nresults, niter = 1000, 1000
    for x in X:
        # We call the fortran function, passing the appropriate Python function.
        results = _chaos.iterate_limit(func(x), 0.5, niter, nresults)
        pl.plot([x]*len(results), results, 'k,')

if __name__ == '__main__':
    pl.figure()
    driver(logistic, 0.0, 4.0)
    pl.xlabel('r')
    pl.ylabel('X limit')
    pl.title('Logistic Map')
    pl.figure()
    driver(sine, 0.0, 1.0)
    pl.xlabel('r')
    pl.ylabel('X limit')
    pl.title('Sine Map')
    pl.show()
