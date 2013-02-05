# Using Python and FORTRAN with F2py 

**Presented by Katy Huff**

**Based on Lecture Material by The Hacker Within**

Motivation: Why use `f2py`?
=============================

Fortran has been around for a while, and it forms the core of much
scientific software. There's much legacy FORTRAN 77 code out there that
works well, is widely tested and optimized. Being FORTRAN 77, no one
wants to touch it unless absolutely necessary.

Python has many advantages over Fortran, and the few areas it is lacking
(primarily execution speed), Python and Fortran are complementary. Where
Python is streamlined and easy to read, Fortran is klunky and prone to
spaghetti code. Python is dynamic, Fortran types everything, even your
spelling mistakes if you're not careful (implicit typing). Python is
young and spry, Fortran is old, pedigreed and barnacled. Python will
nicely tell you exactly where something is going wrong with an exception
and backtrace; often the best you can hope for from Fortran is a seg
fault or bus error, but data corruptions are common and easy to do.

For all its shortcomings, Fortran code -- especially scientific code
with many arrays, etc -- is usually much faster than the Python
equivalent. It would be wonderful to harness the vast libraries of
tested, fast and well-used Fortran code with a sleek Python wrapper.

And that's where `f2py` comes in.

`f2py` is a Fortran to Python interface generator that allows you to
take Fortran code, generate a Python interface and compile it all
together into an extension module. The front-end is Python, the heavy
lifting is done by Fortran. Best of both worlds.

First Step : Getting Help
=========================

Always, we first learn how to get more help. To get some clues about
`f2py`, type f2py in the terminal and you will see hints about usage.

Among the hints about usage, you will find that there is no ordinary
help flag, though there are help flags concerning fortran compilers and
linking. This is good to know, but not super helpful.

Really, help for f2py can mostly just be found on the web :

-   The User's Guide is your friend : http://cens.ioc.ee/projects/f2py2e/usersguide/
-   Google is also your friend. You'll find that f2py might throw errors
    you don't understand. Google them and you'll find the StackOverflow
    forums are full of people with the same error.

Example 0: Hello World
======================

As usual, we want to start with a very simple example provided by the
user manual.

```fortran
C FILE: FIB1.F
      SUBROUTINE FIB(A,N)
C
C     CALCULATE FIRST N FIBONACCI NUMBERS
C
      INTEGER N
      REAL*8 A(N)
      DO I=1,N
         IF (I.EQ.1) THEN
            A(I) = 0.0D0
         ELSEIF (I.EQ.2) THEN
            A(I) = 1.0D0
         ELSE 
            A(I) = A(I-1) + A(I-2)
         ENDIF
      ENDDO
      END
C END FILE FIB1.F
```

The purpose of this simple file is to fill the array you provide with a
fibonacci series.

In order to Pythonize this code, try:

    $ f2py -c -m fib1 fib1.f

The configuration of my machine requires that I specify the compiler I want
to use, so the command that I'll call is : `f2py -c -m --fcompiler=gnu95
fib1 fib1.f` Once we have run this command, a shared object file has been
created by f2py.

```python
  import fib1
```
Interestingly, if we don't yet know how to use fib1 or the fib module
within it , we can view the docstrings created by f2py.

```python
  print fib1.__doc__
  print fib1.fib.__doc__
```

So, now we know that in order to use the fib code, we need to provide a
numpy (Numeric) array to fill with fibonacci numbers.

```python
  import numpy as np
  a=np.zeros(10,'d') 
  fib1.fib(a) 
  print a 
```

Example 1: passing scalar arguments
===================================

Let's try a more interesting example. Suppose we have an implicitly
typed FORTRAN 77 function that takes a number of scalar arguments. This
might be a subroutine in a legacy FORTRAN 77 code, for example.

```fortran
      subroutine scalar_args(int_in, real_in, int_inout, real_inout,
     \ int_out, real_out)
C This doesn't do anything interesting, just for illustration.
      int_inout = int_in
      real_inout = real_in
      int_out = int_inout
      real_out = real_inout
        
      end subroutine scalar_args
```

It is easy to wrap this subroutine with `f2py`. Here's how.

First, it is necessary to tell `f2py` the intent of each subroutine
argument. `f2py` provides multiple ways to specify how to generate the
interface -- the easiest is to put `f2py`-specific comments right in
the FORTRAN code.

[scalar_args.f](https://raw.github.com/thehackerwithin/PyTrieste/master/f2py/scalar_args.f)

```fortran
      subroutine scalar_args(int_in, real_in, int_inout, real_inout,
     \ int_out, real_out)
C Here are the f2py-specific comments.
Cf2py intent(in) :: int_in, real_in }}
Cf2py intent(inout) :: int_inout, real_inout
Cf2py intent(out) :: int_out, real_out

      int_inout = int_in
      real_inout = real_in
      int_out = int_inout
      real_out = real_inout
        
      end subroutine scalar_args
```

You'll notice that the intent specifications are very similar to Fortran
90-style intent statements. The `f2py` specific comments start with
`Cf2py` for FORTRAN 77 code, and `!f2py` for Fortran 9x code.

These intent specifications are necessary for `f2py` to generate the
correct interface. If you're writing Fortran 9x code with intent
specifiers already in place, `f2py` will take care of this for you.

To create the extension module, we invoke `f2py` from the command
line. On UNIX/Linux, assuming the above subroutine is in a source file
'scalar_args.f':

    $ f2py -c -m _scalar_args scalar_args.f

The '-c' switch tells `f2py` to compile an extension module, and the
'-m _scalar_args' specifies the name of the extension module. The
fortran source files follow (in this case just one file).

If everything is setup correctly, the above command will compile the
fortran sources into an extension module named '_scalar_args.so' (the
extension will be different for Mac OS X or Windows) located in the
current directory.

We can test this module from python with a python source file named
'pass_args.py':

[pass_args.py](https://raw.github.com/thehackerwithin/PyTrieste/master/f2py/pass_args.py)

```python
  # pass_args.py
  import numpy as np
  import _scalar_args
  
  print _scalar_args.scalar_args.__doc__
  
  # these are simple python scalars.
  int_in = 1.0
  real_in = 10.0
  
  # since these are intent(inout) variables, these must be arrays
  int_inout = np.zeros((1,), dtype = np.int32)
  real_inout = np.zeros((1,), dtype = np.float32)
  
  # all intent(out) variables are returned in a tuple, so they aren't passed as
  # arguments.
  
  int_out, real_out = _scalar_args.scalar_args(int_in, real_in, int_inout, real_inout)
  
  for name in ('int_inout', 'real_inout', 'int_out', 'real_out'):
      print '%s == %s' % (name, locals()[name])

```

Running the above python script should yield the following output:

    scalar_args - Function signature:
      int_out,real_out = scalar_args(int_in,real_in,int_inout,real_inout)
    Required arguments:
      int_in : input int
      real_in : input float
      int_inout : in/output rank-0 array(int,'i')
      real_inout : in/output rank-0 array(float,'f')
    Return objects:
      int_out : int
      real_out : float

    int=inout == [1]
    real_inout == [ 10.]
    int_out == 1
    real_out == 10.0

One nice feature of `f2py` is that it generates informative docstrings
for the wrapped fortran subroutines & functions. In this case, it tells
us that the subroutine 'scalar_args' has a function signature that
takes 4 inputs and returns a 2-tuple. The first 2 inputs are an int and
a float, respectively. These are the 'intent(in)' variables.

The remaining inputs are 'in/output rank-0 array' objects -- these are
simply numpy arrays with a single element (a rank-0 object). These are
necessary since the fortran objects are intent(inout), and there must be
a place to put the output value.

The intent(out) arguments are converted by `f2py` into a return
2-tuple, and are returned by the wrapper function. This is the case for
any Fortran procedure argument that has an intent(out) attribute.

The remainder of the output shows that the subroutine is behaving
correctly.

Let's move on to passing arrays between Python and Fortran.

Example 2: passing array arguments
==================================

Here's the source of a FORTRAN 77 subroutine that takes array arguments:

```fortran

      subroutine array_args(nx, ny, int_arr_in,
     \ cmplx_arr_inout, 
     \ real_arr_out)

          integer nx, ny
          integer int_arr_in(nx, ny)
          complex cmplx_arr_inout(nx, ny)
          real real_arr_out(nx, ny)

          integer i, j

          do j = 1, ny
              do i = 1, nx
                  cmplx_arr_inout(i,j) = cmplx(int_arr_in(i,j),
     \                   int_arr_in(i,j))
                  real_arr_out(i,j) = real(int_arr_in(i,j))
              enddo
          enddo

      end subroutine array_args
```

Nothing special. This contrived example is designed to be similar to
FORTRAN 77 legacy code that has array arguments, with the array extents
passed in explicitly. You should note that in the loop, the arrays are
iterated through in column-major order (i.e. the first index varies the
fastest). This is known in NumPy & `f2py` parlance as 'fortran order'.
We'll have to keep this in mind when passing multi-dimensional arrays
between Python and Fortran, since Python uses row-major ordering, known
as 'C order'. For 2-dimensional arrays, the orderings are the transpose
of each other, and to index the same element, the array indices need to
be reversed.

Let's add in the `f2py` comments to specify the intent of the
arguments:

[array_args.f](https://raw.github.com/thehackerwithin/PyTrieste/master/f2py/array_args.f)

```fortran
      subroutine array_args(nx, ny, int_arr_in,
     \                      cmplx_arr_inout, 
     \                      real_arr_out)

          integer nx, ny
          integer int_arr_in(nx, ny)
          complex cmplx_arr_inout(nx, ny)
          real real_arr_out(nx, ny)

Cf2py intent(in) nx, ny
Cf2py intent(in) int_arr_in
Cf2py intent(inout) cmplx_arr_inout
Cf2py intent(out) real_arr_out

C ... body of subroutine unchanged ...

      end subroutine array_args

```

As you'd expect. We invoke `f2py` from the commandline:

    $ f2py -c -m _array_args array_args.f

Here's a test script similar to the one we saw before:

[pass_array_args.py](https://raw.github.com/thehackerwithin/PyTrieste/master/f2py/pass_array_args.py)

```python
  # pass_array_args.py
  import numpy as np
  import _array_args
  
  print _array_args.array_args.__doc__
  
  # int_arr is a 10 X 10 array filled with consecutive integers.
  # It is in 'fortran' order.
  int_arr = np.asfortranarray(np.arange(100, dtype = 'i').reshape(10,10))
  
  # cplx_arr is a 10 X 10 complex array filled with zeros.
  # It is in 'fortran' order.
  cplx_arr = np.asfortranarray(np.zeros((10,10), dtype = 'F'))
  
  # We invoke the wrapped fortran subroutine.
  real_arr = _array_args.array_args(int_arr, cplx_arr)
  
  # Here are the results.
  print "int_arr  = %s" %  int_arr
  print "real_arr = %s" % real_arr
  print "cplx_arr = %s" % cplx_arr

```

One thing to note here: the `int_arr` and `cplx_arr` are declared
as fortran arrays, (`np.asfortranarray(...)`) since that's what we
want in this case. Their memory layout is fortran contiguous, and the
fortran subroutine won't have any complaints.

The docstring for the wrapped subroutine is again very helpful:

    array_args - Function signature:
      real_arr_out = array_args(int_arr_in,cmplx_arr_inout,[nx,ny])
    Required arguments:
      int_arr_in : input rank-2 array('i') with bounds (nx,ny)
      cmplx_arr_inout : in/output rank-2 array('F') with bounds (nx,ny)
    Optional arguments:
      nx := shape(int_arr_in,0) input int
      ny := shape(int_arr_in,1) input int
    Return objects:
      real_arr_out : rank-2 array('f') with bounds (nx,ny)

The docstring tells us the subroutine takes 2 arguments, the first a
rank-2 integer array, the second a rank-2 complex array (that's the
`array('F')` part). It is unnecessary to pass in the array extents
explicitly, since the extents can be queried `f2py` from the numpy
arrays themselves.

It also tells us the shape and type of the return array.

The script output gives us the following:

    int_arr  == [[ 0  1  2  3  4  5  6  7  8  9]
     [10 11 12 13 14 15 16 17 18 19]
     [20 21 22 23 24 25 26 27 28 29]
     [30 31 32 33 34 35 36 37 38 39]
     [40 41 42 43 44 45 46 47 48 49]
     [50 51 52 53 54 55 56 57 58 59]
     [60 61 62 63 64 65 66 67 68 69]
     [70 71 72 73 74 75 76 77 78 79]
     [80 81 82 83 84 85 86 87 88 89]
     [90 91 92 93 94 95 96 97 98 99]]
    real_arr == [[  0.   1.   2.   3.   4.   5.   6.   7.   8.   9.]
     [ 10.  11.  12.  13.  14.  15.  16.  17.  18.  19.]
     [ 20.  21.  22.  23.  24.  25.  26.  27.  28.  29.]
     [ 30.  31.  32.  33.  34.  35.  36.  37.  38.  39.]
     [ 40.  41.  42.  43.  44.  45.  46.  47.  48.  49.]
     [ 50.  51.  52.  53.  54.  55.  56.  57.  58.  59.]
     [ 60.  61.  62.  63.  64.  65.  66.  67.  68.  69.]
     [ 70.  71.  72.  73.  74.  75.  76.  77.  78.  79.]
     [ 80.  81.  82.  83.  84.  85.  86.  87.  88.  89.]
     [ 90.  91.  92.  93.  94.  95.  96.  97.  98.  99.]]
    cplx_arr == [[  0. +0.j   1. +1.j   2. +2.j   3. +3.j   4. +4.j   5. +5.j   6. +6.j
        7. +7.j   8. +8.j   9. +9.j]
     [ 10.+10.j  11.+11.j  12.+12.j  13.+13.j  14.+14.j  15.+15.j  16.+16.j
       17.+17.j  18.+18.j  19.+19.j]
     [ 20.+20.j  21.+21.j  22.+22.j  23.+23.j  24.+24.j  25.+25.j  26.+26.j
       27.+27.j  28.+28.j  29.+29.j]
     [ 30.+30.j  31.+31.j  32.+32.j  33.+33.j  34.+34.j  35.+35.j  36.+36.j
       37.+37.j  38.+38.j  39.+39.j]
     [ 40.+40.j  41.+41.j  42.+42.j  43.+43.j  44.+44.j  45.+45.j  46.+46.j
       47.+47.j  48.+48.j  49.+49.j]
     [ 50.+50.j  51.+51.j  52.+52.j  53.+53.j  54.+54.j  55.+55.j  56.+56.j
       57.+57.j  58.+58.j  59.+59.j]
     [ 60.+60.j  61.+61.j  62.+62.j  63.+63.j  64.+64.j  65.+65.j  66.+66.j
       67.+67.j  68.+68.j  69.+69.j]
     [ 70.+70.j  71.+71.j  72.+72.j  73.+73.j  74.+74.j  75.+75.j  76.+76.j
       77.+77.j  78.+78.j  79.+79.j]
     [ 80.+80.j  81.+81.j  82.+82.j  83.+83.j  84.+84.j  85.+85.j  86.+86.j
       87.+87.j  88.+88.j  89.+89.j]
     [ 90.+90.j  91.+91.j  92.+92.j  93.+93.j  94.+94.j  95.+95.j  96.+96.j
       97.+97.j  98.+98.j  99.+99.j]]

What if we had not declared the `int_arr` as fortran contiguous?
Let's see.

First, let's turn-on array-copying output in the fortran module. This
requires us to recompile the module with a commandline flag.

    $ f2py -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m _array_args array_args.f

Let's change the pass_array_args.py file thusly:

```python
  # int_arr is a 10 X 10 array filled with consecutive integers.
  # Now it is in 'C' order.
  int_arr = np.arange(100, dtype = 'i').reshape(10,10)
```

When running the script, you will notice an extra line in the output:

```
  copied an array: size = 100, elsize = 4
  int_arr  == [[ 0  1  2  3  4  5  6  7  8  9]
  ...
  real_arr == [[  0.   1.   2.   3.   4.   5.   6.   7.   8.   9.]
  ...
  cplx_arr == [[  0. +0.j   1. +1.j   2. +2.j   3. +3.j   4. +4.j   5. +5.j   6. +6.j
  ...

```

The `-DF2PY_REPORT_ON_ARRAY_COPY=1` switch caused `f2py` to
report that it copied an array (int_arr) on input, since it received a
'C' order array as an argument. To avoid this array copy, it is
necessary to declare the arrays as fortran contiguous, with the
`np.asfortranarray` function.

Example 3: using .pyf files and Python callbacks
================================================

The above 2 examples, while simple, cover a large chunk of calling
Fortran from Python with `f2py`. It is possible to call Python from
Fortran, using callbacks.

As a more interesting example, we'll plot the [logistic
map](http://en.wikipedia.org/wiki/Logistic_map), a classic plot exhibiting
self-similarity and period-doubling yielding chaos and fractal structure.
The logistic map is a fascinating system that shows how very simple
nonlinear systems have nearly unlimited richness. It can be used as a very
simple model of year-to-year populations that are limited by resources or
subject to predator-prey dynamics (I'm a plasma physicist, not an
ecologist, so don't harangue me over the details!).

Let's say we have a Fortran subroutine that calculates the equilibrium
points for an iteratively applied function. It takes a function as an
argument, applies the function iteratively `num_iters` times, and
puts the next `n` results of the function in an array.

[chaos.f](https://raw.github.com/thehackerwithin/PyTrieste/master/f2py/chaos.f)

```fortran
      subroutine iterate_limit(func, x0, num_iters, results, n)
          external func
          double precision func
          double precision x0
          integer num_iters, n
          double precision results(n)

          integer i

          do i = 1, num_iters
              x0 = func(x0)
          enddo

          do i = 1, n
              results(i) = x0
              x0 = func(x0)
          enddo

      end subroutine iterate_limit
```

The above is saved in a file `chaos.f`.

This time, rather than put `Cf2py` comments in the Fortran source,
we'll instead generate an interface file.

Call `f2py` thusly:

    $ f2py -h _chaos.pyf chaos.f

This command instructs `f2py` to extract the necessary information
from the fortran source and create an interface file `_chaos.pyf`
that we'll edit accordingly.

Here's the output:

```fortran

!    -*- f90 -*-
! Note: the context of this file is case sensitive.

python module iterate_limit__user__routines 
    interface iterate_limit_user_interface 
        function func(x0) result (x0) ! in :_chaos:chaos.f:iterate_limit:unknown_interface
            double precision :: x0
        end function func
    end interface iterate_limit_user_interface
end python module iterate_limit__user__routines
python module _chaos ! in 
    interface  ! in :_chaos
        subroutine iterate_limit(func,x0,num_iters,results,n) ! in :_chaos:chaos.f
            use iterate_limit__user__routines
            external func
            double precision :: x0
            integer :: num_iters
            double precision dimension(n) :: results
            integer optional,check(len(results)> = n),depend(results) :: n = len(results)
        end subroutine iterate_limit
    end interface 
end python module _chaos

! This file was auto-generated with f2py (version:2).
! See http://cens.ioc.ee/projects/f2py2e/

```

All that remains in this instance is to add in intent specifications to the
interface file, and to remove the line specifying the `n` argument. Here
are the changed lines
[chaos.pyf](https://raw.github.com/thehackerwithin/PyTrieste/master/f2py/chaos.pyf):

```python
            double precision, intent(inout) :: x0
            integer, intent(in) :: num_iters
            integer, intent(in) :: n
            double precision dimension(n), intent(out) :: results
```

Now, we invoke `f2py` a bit differently, to use the interface file.

    $ f2py -c -m _chaos _chaos.pyf chaos.f

Here is the driver script in Python
[chaos.py](https://raw.github.com/thehackerwithin/PyTrieste/master/f2py/chaos.py):

```python
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
```

Running the above script yields [attachment:logistic-map.png] and
[attachment:sine-map.png] .

[[Image(logistic-map.png, 50%, center, top)]]

[[Image(sine-map.png, 50%, center, top)]]

The significance to note here is that we are able to pass an arbitrary
Python function (provided it has the right signature!) to Fortran code,
the Fortran code calls the Python function and does something useful
with it. We can easily change which function is passed from Python, thus
achieving a greater degree of flexibility using Python with Fortran.

Conclusions & External links
============================

There's much more to `f2py` than presented here -- here are some
useful links.

 - The `f2py` documentation -- http://cens.ioc.ee/projects/f2py2e/
 - `f2py` is distributed as part of NumPy --
   [http://numpy.scipy.org](http://numpy.scipy.org)/
 - And `f2py` is used to generate much of the wrappers for SciPy --
   [http://www.scipy.org](http://www.scipy.org)/

`f2py` supports some Fortran 9x specific features, and it is possible
to wrap module procedures with `f2py`. Derived types are not
supported, however. Neither are assumed-shape arrays. In short, `f2py`
excels at wrapping FORTRAN 77 code and supports everything any sane
person would want to do with Python and FORTRAN 77.
