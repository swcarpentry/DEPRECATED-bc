Python Testing Cheat Sheet
==========================

Why testing?
------------

1. Helps you to think about expected behavior, especially boundary cases,
2. documents expected behavior,
3. confidence recent changes didn't break anything that worked before,
4. confidence code is correct.


Defensive programming
---------------------

Using an assertion to ensure input is acceptable:

    def some_function(x):
        assert x >= 0
        # ... continue safe in knowledge that x > 0

Adding an explanatory message to the assertion:

        assert x >= 0, "Function not defined for negative x."

Alternatively, raise an exception to indicate what the problem is:

    def some_function(x):
        if x < 0:
            raise TypeError, "Function not defined for negative x."
        return 0


Unit testing with Nose
----------------------

To run tests, at the shell prompt, type

    nosetests

By default, Nose will

* look for test functions that have a name starting with `test`
* look for them in files with names starting with `test`
* look for such files in the current working directory, and in subdirectories with names starting with `test`

There are some additional rules, and you can configure your own, but this should be enough to get started.

### A simple test

    from nose.tools import assert_equal

    from mystatscode import mean

    def test_single_value():
        observed = mean([1.0])
        expected = 1.0
        assert_equal(observed, expected)

### Other assertions

Nose provides a range of assertions that can be used when a test is not just checking a simple equality, e.g.

    from nose.tools import assert_items_equal

    from mycode import find_factors

    def test_6():
        observed = find_factors(6)
        expected = [2, 3]
        assert_items_equal(observed, expected) # order of factors is not guaranteed

* assertTrue, assertFalse
* assertIn, assertNotIn
* assertIs, assertIsNot
* assertRaises
* (what else?)

### Floating point tests

When comparing floating-point numbers for equality, allow some tolerance for small differences due to
the way values are represented and rounded.
* assertGreater, assertLess

    from nose.tools import assert_almost_equal

    from mycode import hypotenuse

    def test_hypotenuse_345():
        observed = hypotenuse(3.0, 4.0)
        expected = 5.0
        assert_almost_equal(observed, expected)

### Testing exceptions

Testing that a method raises the appropriate exception when the input is invalid:

    from nose.tools import raises

    from mystatscode import mean

    @raises(TypeError)
    def test_not_a_list():
        observed = mean(1)

### Fixtures

A *fixture* is what the test function uses as input, e.g. values, objects and arrays.

To set up a fixture once before any tests are run, define a method called `setup` in the same files
as the test functions. This can assign values to global variables for use in the test functions.

    long_list = None

    def setup():
        long_list = [0]
        # append more values to long_list...

If the global variables assigned in `setup` might be modified by some of the test functions, the set-up
step must be executed once before each test function is called:

    from nose.tools import with_setup

    from mycode import mean, clear

    long_list = None

    def setup_each():
        long_list = [0]
        # append more values to long_list...

    @with_setup(setup_each)
    def test_mean_long_list():
        observed = mean(long_list)
        expected = 0.0
        assert_equal(observed, expected)

    @with_setup(setup_each)
    def test_clear_long_list():
        clear(long_list)
	assert_equal(len(long_list), 0)



Test-driven deveopment
----------------------

***Red.*** Write test function that checks one new functionality you want to add to your code. -- tests have to fail.

***Green.*** Write minimal code that implements desired features until all tests pass.

***Refactor.*** Improve code wrt. readability and speed. Constantly check that tests still pass.

***Commit.*** Commit working code to version control.

Repeat.


General advice
--------------

* Perfect test-case coverage is impossible.
* Try to test distinct functionalities.
* If you find a bug yet undiscovered by previous test, make it a new test case.



