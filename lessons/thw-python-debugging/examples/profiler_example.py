#!/usr/bin/env python
"""
A script that shows off profiling via list iteration schemes
"""

from random import random


def dynamic_array(size=1000000):
    """Fills an array that is sized dynamically."""
    dynamic = []
    for i in range(size):
        dynamic.append(random() * i)
    return dynamic


def static_array(size=1000000):
    """Fills an array that is pre-allocated."""
    static = [None] * size
    for i in range(size):
        static[i] = random() * i
    return static


def comprehension_array(size=1000000):
    """Fills an array that is handled by Python via list comprehension."""
    return [random() * i for i in range(size)]

if __name__ == "__main__":
    import sys

    # Allow the user to input filled array sizes
    size = 1000000
    for i, val in enumerate(sys.argv):
        if val == "--size":
            size = int(sys.argv[i + 1])

    # Allow the user to specify the method of array filling
    if "--dynamic" in sys.argv:
        dynamic_array(size)
    if "--static" in sys.argv:
        static_array(size)
    if "--comprehension" in sys.argv:
        comprehension_array(size)
