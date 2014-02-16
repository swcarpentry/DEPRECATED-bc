---
layout: lesson
root: ../..
title: Unit Testing
level: novice
---
Most people don't enjoy writing tests,
so if we want them to actually do it,
it must be easy to:

<ul>
<li>add or change tests,</li>
<li>understand the tests that have already been written,</li>
<li>run those tests, and</li>
<li>understand those tests' results.</li>
</ul>

Test results must also be reliable.
If a testing tool says that code is working when it's not,
or reports problems when there actually aren't any,
people will lose faith in it and stop using it.

The simplest kind of test is a [unit test](../gloss.html#unit-test)
that checks the behavior of one component of a program.
As an example,
suppose we're testing a function called `rectangle_area`
that returns the area of an `(x0, y0, x1, y1)` rectangle.
We'll start by testing our code directly using `assert`.
Here,
we call the function three times with different arguments,
checking that the right value is returned each time.

~~~
from rectangle import rectangle_area

assert rectangle_area([0, 0, 1, 1]) == 1.0
assert rectangle_area([1, 1, 4, 4]) == 9.0
assert rectangle_area([0, 1, 4, 7]) == 24.0

<span class="err">---------------------------------------------------------------------------
AssertionError                            Traceback (most recent call last)

<ipython-input-16-ebf7f5f1c120> in <module>()
3 assert rectangle_area([0, 0, 1, 1]) == 1.0
4 assert rectangle_area([1, 1, 4, 4]) == 9.0
----> 5 assert rectangle_area([0, 1, 4, 7]) == 24.0

AssertionError:
~~~

This result is used,
in the sense that we know something's wrong,
but look closely at what happens if we run the tests in a different order:

~~~
assert rectangle_area([0, 1, 4, 7]) == 24.0
assert rectangle_area([1, 1, 4, 4]) == 9.0
assert rectangle_area([0, 0, 1, 1]) == 1.0

<span class="err">---------------------------------------------------------------------------
AssertionError                            Traceback (most recent call last)

<ipython-input-17-548f3f32c981> in <module>()
----> 1 assert rectangle_area([0, 1, 4, 7]) == 24.0
2 assert rectangle_area([1, 1, 4, 4]) == 9.0
3 assert rectangle_area([0, 0, 1, 1]) == 1.0

AssertionError:
~~~

Python halts at the first failed assertion,
so the second and third tests aren't run at all.
It would be more helpful if we could get data from all of our tests every time
they're run,
since the more information we have,
the faster we're likely to be able to track down bugs.
It would also be helpful to have some kind of summary report:
if our [test suite](../gloss.html#test-suite) includes thirty or forty tests
(as it well might for a complex function or library that's widely used),
we'd like to know how many passed or failed.

Here's a different approach.
First, let's put each test in a function with a meaningful name:

~~~
def test_unit_square():
    assert rectangle_area([0, 0, 1, 1]) == 1.0

def test_large_square():
    assert rectangle_area([1, 1, 4, 4]) == 9.0

def test_actual_rectangle():
    assert rectangle_area([0, 1, 4, 7]) == 24.0
~~~

Next,
import a library called `ears`
and ask it to run our tests for us:

~~~
import ears
ears.run()

..f
2 pass, 1 fail, 0 error
----------------------------------------
fail: test_actual_rectangle
Traceback (most recent call last):
  File "ears.py", line 45, in run
test()
  File "<ipython-input-18-643689ad0a0f>", line 8, in test_actual_rectangle
assert rectangle_area([0, 1, 4, 7]) == 24.0
AssertionError
~~~

`ears.run` looks in the calling program
for functions whose names start with the letters `'test_'`
and runs each one.
If the function complete without an assertion being triggered,
we count the test as a [success](../gloss.html#test-success).
If an assertion fails,
we count the test as a [failure](../gloss.html#test-failure),
but if any other exception occurs,
we count it as an [error](../gloss.html#test_error)
because the odds are that the test itself is broken.

`ears` is an [xUnit testing library](../gloss.html#xunit).
The name "xUnit" comes from the fact that
many of them are imitations of a Java testing library called JUnit.
The [Wikipedia page](http://en.wikipedia.org/wiki/List_of_unit_testing_frameworks) on the subject subject
lists dozens of similar frameworks in almost as many languages,
all of which have a similar structure:
each test is a single function that follows some naming convention
(e.g., starts with `'test_'`),
and the framework runs them in some order
and reports how many passed, failed, or were broken.
